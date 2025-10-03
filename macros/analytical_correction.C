// analytical_correction.C — purely analytical Λ→pπ acceptance and self-made plots
// Usage: root -l -q analytical_correction.C

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TFile.h>
#include <TSystem.h>

#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>

namespace ana {

// ------------------------- Physical constants (GeV, cm) -------------------------
constexpr double mL   = 1.115683;  // Λ mass
constexpr double mp   = 0.938272;  // proton mass
constexpr double mpi  = 0.139570;  // charged pion mass
constexpr double ctau = 7.89;      // cτ_Λ (cm)
constexpr double BpPi = 0.639;     // Br(Λ→pπ−)
constexpr double alphaLambda = 0.732; // PDG decay asymmetry α_Λ

// ------------------------- Configurable knobs -----------------------------------
struct AcceptCfg {
  double Lmin_cm   = 1.5;   // resolvable PV–DV separation
  double Lmax_cm   = 200.;  // available in-volume flight (set large to ignore)
  double pthr_p    = 0.25;  // proton p-threshold [GeV/c]
  double pthr_pi   = 0.10;  // pion  p-threshold [GeV/c]
  double PLambda   = 0.0;   // longitudinal polarisation (0 = unpolarised)
  double alpha     = alphaLambda;
};
static AcceptCfg gCfg{};

// ------------------------- Helper kinematics (Eqs. 6.3–6.6) ---------------------
inline double kallen(double a, double b, double c) {
  return a*a + b*b + c*c - 2.0*(a*b + a*c + b*c);
}
inline double beta_from_p(double p) {
  const double E = std::sqrt(p*p + mL*mL);
  return (E > 0) ? p / E : 0.0;
}
inline double gamma_from_p(double p) {
  return std::sqrt(1.0 + (p*p)/(mL*mL));
}

// Alength: e^{−Lmin/λ} − e^{−Lmax/λ}, λ = βγ cτ_Λ  (Eq. 6.2)
inline double Alength(double beta_gamma, double Lmin, double Lmax) {
  const double lambda = std::max(1e-12, beta_gamma * ctau);
  const double t1 = std::exp(-Lmin / lambda);
  const double t2 = (Lmax > 0 ? std::exp(-Lmax / lambda) : 0.0);
  return std::clamp(t1 - t2, 0.0, 1.0);
}

// Akin: two-body thresholds + optional polarisation (Eqs. 6.4–6.8)
inline double Akin(double p, double pthr_p, double pthr_pi, double P, double alpha) {
  const double beta  = beta_from_p(p);
  const double gamma = gamma_from_p(p);
  if (beta <= 0.0) return 0.0;

  const double lam = kallen(mL*mL, mp*mp, mpi*mpi);
  if (lam <= 0.0) return 0.0;

  const double pst   = 0.5 * std::sqrt(lam) / mL;     // p* in Λ rest frame (Eq. 6.3)
  const double Epst  = std::sqrt(mp*mp  + pst*pst);
  const double Epist = std::sqrt(mpi*mpi + pst*pst);

  // Implement p-thresholds as energy thresholds (Eq. 6.5)
  const double Ethr_p  = std::sqrt(mp*mp  + pthr_p*pthr_p);
  const double Ethr_pi = std::sqrt(mpi*mpi + pthr_pi*pthr_pi);

  const double denom = beta * pst; // linear bounds in cosθ* (Eq. 6.6)
  if (denom <= 0.0) return 0.0;

  const double cmin_p  = (Ethr_p/gamma  - Epst) / denom; // x ≥ cmin_p
  const double cmax_pi = (Epist - Ethr_pi/gamma) / denom; // x ≤ cmax_pi

  double l = std::max(-1.0, cmin_p);
  double u = std::min( 1.0, cmax_pi);
  if (l >= u) return 0.0;

  // Unpolarised fraction plus polarisation correction (Eqs. 6.7–6.8)
  const double Aiso = 0.5*(u - l);
  const double Apol = 0.25*alpha*P*(u*u - l*l);
  return std::clamp(Aiso + Apol, 0.0, 1.0);
}

// APS: B × Alength × Akin (Eq. 6.1)
inline double APS_from_p(double p) {
  const double beta_gamma = p / mL; // βγ = p/m
  const double Alen = Alength(beta_gamma, gCfg.Lmin_cm, gCfg.Lmax_cm);
  const double A2   = Akin(p, gCfg.pthr_p, gCfg.pthr_pi, gCfg.PLambda, gCfg.alpha);
  const double aps  = BpPi * Alen * A2;
  return std::clamp(aps, 0.0, 1.0);
}
inline double APS_from_bg(double beta_gamma) { // convenience: APS vs βγ
  const double p = beta_gamma * mL;
  return APS_from_p(p);
}

} // namespace ana

// ------------------------- Plot builders (self-made) ------------------------------
void draw_APS_vs_bg(const std::string& outdir) {
  const int N = 600;
  std::vector<double> x(N), y_nom(N), y_loL(N), y_hiL(N), y_noThr(N), y_polp(N), y_polm(N);

  // Variants: default Lmin, low Lmin, high Lmin; no thresholds; PΛ = ±0.4
  ana::AcceptCfg cfg_nom = ana::gCfg;
  ana::AcceptCfg cfg_loL = ana::gCfg; cfg_loL.Lmin_cm = 0.5;
  ana::AcceptCfg cfg_hiL = ana::gCfg; cfg_hiL.Lmin_cm = 2.5;
  ana::AcceptCfg cfg_thr = ana::gCfg; cfg_thr.pthr_p = 0.0; cfg_thr.pthr_pi = 0.0;
  ana::AcceptCfg cfg_pP  = ana::gCfg; cfg_pP.PLambda = +0.4;
  ana::AcceptCfg cfg_mP  = ana::gCfg; cfg_mP.PLambda = -0.4;

  auto eval_APS = [](double bg, const ana::AcceptCfg& C)->double {
    const double p = bg * ana::mL;
    const double Alen = ana::Alength(bg, C.Lmin_cm, C.Lmax_cm);
    const double A2   = ana::Akin(p, C.pthr_p, C.pthr_pi, C.PLambda, C.alpha);
    return std::clamp(ana::BpPi * Alen * A2, 0.0, 1.0);
  };

  const double bg_min = 0.0, bg_max = 10.0;
  for (int i = 0; i < N; ++i) {
    const double bg = bg_min + (bg_max - bg_min) * i / (N - 1);
    x[i]      = bg;
    y_nom[i]  = eval_APS(bg, cfg_nom);
    y_loL[i]  = eval_APS(bg, cfg_loL);
    y_hiL[i]  = eval_APS(bg, cfg_hiL);
    y_noThr[i]= eval_APS(bg, cfg_thr);
    y_polp[i] = eval_APS(bg, cfg_pP);
    y_polm[i] = eval_APS(bg, cfg_mP);
  }

  TCanvas c("c_APS","APS vs beta*gamma",900,700);
  gStyle->SetOptStat(0);
  auto g_nom   = new TGraph(N, x.data(), y_nom.data());
  auto g_loL   = new TGraph(N, x.data(), y_loL.data());
  auto g_hiL   = new TGraph(N, x.data(), y_hiL.data());
  auto g_thr   = new TGraph(N, x.data(), y_noThr.data());
  auto g_polp  = new TGraph(N, x.data(), y_polp.data());
  auto g_polm  = new TGraph(N, x.data(), y_polm.data());

  g_nom ->SetLineColor(kBlue+1);  g_nom ->SetLineWidth(3);
  g_loL ->SetLineColor(kGreen+2); g_loL ->SetLineStyle(2); g_loL ->SetLineWidth(2);
  g_hiL ->SetLineColor(kGreen+2); g_hiL ->SetLineStyle(9); g_hiL ->SetLineWidth(2);
  g_thr ->SetLineColor(kRed+1);   g_thr ->SetLineStyle(3); g_thr ->SetLineWidth(2);
  g_polp->SetLineColor(kMagenta+1); g_polp->SetLineStyle(7); g_polp->SetLineWidth(2);
  g_polm->SetLineColor(kMagenta+1); g_polm->SetLineStyle(7); g_polm->SetLineWidth(2);

  g_nom->SetTitle(";#beta#gamma ( #it{#Lambda} boost );A_{PS}(#beta#gamma)");
  g_nom->GetXaxis()->SetLimits(bg_min, bg_max);
  g_nom->GetYaxis()->SetRangeUser(0.0, 0.8);

  g_nom->Draw("AL");
  g_thr->Draw("L SAME");
  g_loL->Draw("L SAME");
  g_hiL->Draw("L SAME");
  g_polp->Draw("L SAME");
  g_polm->Draw("L SAME");

  TLegend leg(0.50,0.18,0.88,0.40);
  leg.SetBorderSize(0);
  leg.AddEntry(g_nom , Form("Nominal: L_{min}=%.1f cm, p_{thr}^{p}=%.2f GeV, p_{thr}^{#pi}=%.2f GeV, P_{#Lambda}=0",
                            ana::gCfg.Lmin_cm, ana::gCfg.pthr_p, ana::gCfg.pthr_pi), "l");
  leg.AddEntry(g_thr , "No daughter thresholds", "l");
  leg.AddEntry(g_loL , "L_{min}=0.5 cm", "l");
  leg.AddEntry(g_hiL , "L_{min}=2.5 cm", "l");
  leg.AddEntry(g_polp, "P_{#Lambda}=+0.4", "l");
  leg.AddEntry(g_polm, "P_{#Lambda}=-0.4", "l");
  leg.Draw();

  TLatex latex;
  latex.SetNDC(); latex.SetTextSize(0.035);
  latex.DrawLatex(0.13,0.92,"Analytical acceptance for #Lambda #rightarrow p#pi^{-}:  A_{PS}=B #times A_{length} #times A_{kin}");

  c.SaveAs((outdir + "/APS_vs_bg.png").c_str());
  c.SaveAs((outdir + "/APS_vs_bg.pdf").c_str());
}

void draw_Akin_vs_bg(const std::string& outdir) {
  // Show purely the kinematic factor (thresholds & polarisation) versus βγ
  const int N = 600;
  const double bg_min = 0.0, bg_max = 10.0;
  std::vector<double> x(N), y_nom(N), y_noThr(N), y_polp(N), y_polm(N);

  for (int i = 0; i < N; ++i) {
    const double bg = bg_min + (bg_max - bg_min) * i / (N - 1);
    const double p = bg * ana::mL;
    x[i]      = bg;
    y_nom[i]  = ana::Akin(p, ana::gCfg.pthr_p, ana::gCfg.pthr_pi, 0.0, ana::gCfg.alpha);
    y_noThr[i]= ana::Akin(p, 0.0, 0.0, 0.0, ana::gCfg.alpha);
    y_polp[i] = ana::Akin(p, ana::gCfg.pthr_p, ana::gCfg.pthr_pi, +0.4, ana::gCfg.alpha);
    y_polm[i] = ana::Akin(p, ana::gCfg.pthr_p, ana::gCfg.pthr_pi, -0.4, ana::gCfg.alpha);
  }

  TCanvas c("c_Akin","A_kin vs beta*gamma",900,700);
  gStyle->SetOptStat(0);
  auto g_nom   = new TGraph(N, x.data(), y_nom.data());
  auto g_thr   = new TGraph(N, x.data(), y_noThr.data());
  auto g_polp  = new TGraph(N, x.data(), y_polp.data());
  auto g_polm  = new TGraph(N, x.data(), y_polm.data());

  g_nom ->SetLineColor(kBlue+1);  g_nom ->SetLineWidth(3);
  g_thr ->SetLineColor(kRed+1);   g_thr ->SetLineStyle(3); g_thr ->SetLineWidth(2);
  g_polp->SetLineColor(kMagenta+1); g_polp->SetLineStyle(7); g_polp->SetLineWidth(2);
  g_polm->SetLineColor(kMagenta+1); g_polm->SetLineStyle(7); g_polm->SetLineWidth(2);

  g_nom->SetTitle(";#beta#gamma ( #it{#Lambda} boost );A_{kin}");
  g_nom->GetXaxis()->SetLimits(bg_min, bg_max);
  g_nom->GetYaxis()->SetRangeUser(0.0, 1.0);

  g_nom->Draw("AL");
  g_thr->Draw("L SAME");
  g_polp->Draw("L SAME");
  g_polm->Draw("L SAME");

  TLegend leg(0.50,0.18,0.86,0.38);
  leg.SetBorderSize(0);
  leg.AddEntry(g_nom , Form("Nominal thresholds: p_{thr}^{p}=%.2f GeV, p_{thr}^{#pi}=%.2f GeV", ana::gCfg.pthr_p, ana::gCfg.pthr_pi), "l");
  leg.AddEntry(g_thr , "No thresholds", "l");
  leg.AddEntry(g_polp, "Polarised P_{#Lambda}=+0.4", "l");
  leg.AddEntry(g_polm, "Polarised P_{#Lambda}=-0.4", "l");
  leg.Draw();

  TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035);
  latex.DrawLatex(0.13,0.92,"Two-body threshold factor  A_{kin}  (Eqs. 6.4–6.8)");

  c.SaveAs((outdir + "/Akin_vs_bg.png").c_str());
  c.SaveAs((outdir + "/Akin_vs_bg.pdf").c_str());
}

void draw_APS_heatmap_Lmin_vs_bg(const std::string& outdir) {
  // Heatmap: A_PS(βγ, L_min). L_max fixed, thresholds fixed.
  const int Nx = 400, Ny = 120;
  const double bg_min = 0.0, bg_max = 10.0;
  const double Lmin_min = 0.0, Lmin_max = 3.0;

  TH2D h("hAPS",";#beta#gamma ( #it{#Lambda} boost );L_{min} [cm]",
         Nx, bg_min, bg_max, Ny, Lmin_min, Lmin_max);

  // Temporary copy of config to sweep Lmin
  ana::AcceptCfg cfg = ana::gCfg;

  for (int ix = 1; ix <= Nx; ++ix) {
    const double bg = h.GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= Ny; ++iy) {
      cfg.Lmin_cm = h.GetYaxis()->GetBinCenter(iy);
      const double p = bg * ana::mL;
      const double Alen = ana::Alength(bg, cfg.Lmin_cm, cfg.Lmax_cm);
      const double A2   = ana::Akin(p, cfg.pthr_p, cfg.pthr_pi, cfg.PLambda, cfg.alpha);
      const double aps  = std::clamp(ana::BpPi * Alen * A2, 0.0, 1.0);
      h.SetBinContent(ix, iy, aps);
    }
  }

  TCanvas c("c_heat","A_{PS}(#beta#gamma, L_{min})",900,700);
  gStyle->SetOptStat(0);
  h.SetContour(60);
  h.GetZaxis()->SetRangeUser(0.0, 0.7);
  h.Draw("COLZ");
  TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035);
  latex.DrawLatex(0.13,0.92,Form("A_{PS}=B#timesA_{length}#timesA_{kin};  L_{max}=%.0f cm,  p_{thr}^{p}=%.2f GeV,  p_{thr}^{#pi}=%.2f GeV",
                                 ana::gCfg.Lmax_cm, ana::gCfg.pthr_p, ana::gCfg.pthr_pi));
  c.SaveAs((outdir + "/APS_heatmap_Lmin_vs_bg.png").c_str());
  c.SaveAs((outdir + "/APS_heatmap_Lmin_vs_bg.pdf").c_str());
}

// ------------------------- Top-level: analytical_correction() -------------------------------
void analytical_correction() {
  // Configure the analytical acceptance (edit to taste)
  ana::gCfg = {
    .Lmin_cm = 1.5,   // cm
    .Lmax_cm = 200.,  // cm
    .pthr_p  = 0.25,  // GeV/c
    .pthr_pi = 0.10,  // GeV/c
    .PLambda = 0.0,   // unpolarised
    .alpha   = ana::alphaLambda
  };

  const std::string outdir = "plots/analytical";
  gSystem->mkdir(outdir.c_str(), kTRUE);
  gStyle->SetNumberContours(60);

  std::cout << "[analytical] Using: Lmin=" << ana::gCfg.Lmin_cm << " cm, "
            << "Lmax=" << ana::gCfg.Lmax_cm << " cm, "
            << "p_thr^p=" << ana::gCfg.pthr_p << " GeV, "
            << "p_thr^pi=" << ana::gCfg.pthr_pi << " GeV, "
            << "P_Lambda=" << ana::gCfg.PLambda << std::endl;

  draw_APS_vs_bg(outdir);
  draw_Akin_vs_bg(outdir);
  draw_APS_heatmap_Lmin_vs_bg(outdir);

  // Optional: store graphs/heatmap in a ROOT file (comment out if not needed)
  // TFile fout((outdir + "/analytical_acceptance.root").c_str(), "RECREATE");
  // (Recreate the objects here or write from memory if you hold pointers.)
  // fout.Close();
}
