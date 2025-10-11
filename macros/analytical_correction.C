

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "rarexsec/Plotter.hh"

namespace ana {

constexpr double mL = 1.115683;
constexpr double mp = 0.938272;
constexpr double mpi = 0.139570;
constexpr double ctau = 7.89;
constexpr double BpPi = 0.639;
constexpr double alphaLambda = 0.732;

struct AcceptCfg {
    double Lmin_cm = 1.5;
    double Lmax_cm = 200.;
    double pthr_p = 0.25;
    double pthr_pi = 0.10;
    double PLambda = 0.0;
    double alpha = alphaLambda;

    double sigmaE_p = 0.010;
    double sigmaE_pi = 0.010;
};
static AcceptCfg gCfg{};

inline double kallen(double a, double b, double c) {
    return a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c);
}
inline double beta_from_p(double p) {
    const double E = std::sqrt(p * p + mL * mL);
    return (E > 0) ? p / E : 0.0;
}
inline double gamma_from_p(double p) {
    return std::sqrt(1.0 + (p * p) / (mL * mL));
}

inline double Alength(double beta_gamma, double Lmin, double Lmax) {
    const double lambda = std::max(1e-12, beta_gamma * ctau);
    const double t1 = std::exp(-Lmin / lambda);
    const double t2 = (Lmax > 0 ? std::exp(-Lmax / lambda) : 0.0);

    return t1 - t2;
}

inline double Akin(double p, double pthr_p, double pthr_pi, double P, double alpha) {
    const double beta = beta_from_p(p);
    const double gamma = gamma_from_p(p);
    if (beta <= 0.0)
        return 0.0;

    const double lam = kallen(mL * mL, mp * mp, mpi * mpi);
    if (lam <= 0.0)
        return 0.0;

    const double pst = 0.5 * std::sqrt(lam) / mL;
    const double Epst = std::sqrt(mp * mp + pst * pst);
    const double Epist = std::sqrt(mpi * mpi + pst * pst);

    const double Ethr_p = std::sqrt(mp * mp + pthr_p * pthr_p);
    const double Ethr_pi = std::sqrt(mpi * mpi + pthr_pi * pthr_pi);

    const double denom = beta * pst;
    if (denom <= 0.0)
        return 0.0;

    const double cmin_p = (Ethr_p / gamma - Epst) / denom;
    const double cmax_pi = (Epist - Ethr_pi / gamma) / denom;
    const double raw_width = cmax_pi - cmin_p;
    const double tol = 64 * std::numeric_limits<double>::epsilon() *
                       std::max({std::abs(cmin_p), std::abs(cmax_pi), 1.0});
    if (raw_width <= tol)
        return 0.0;

    const double l = std::clamp(cmin_p, -1.0, 1.0);
    const double u = std::clamp(cmax_pi, -1.0, 1.0);
    const double width = std::max(u - l, 0.0);
    if (width <= 0.0)
        return 0.0;

    const double Aiso = 0.5 * width;
    const double Apol = 0.25 * alpha * P * (u * u - l * l);
    return Aiso + Apol;
}

static constexpr double GL64_X[64] = {
    -0.99930504173577217042, -0.99634011677195521983, -0.99101337147674428696, -0.98333625388462597705,
    -0.97332682778991097550, -0.96100879965205376898, -0.94641137485840276522, -0.92956917213193956950,
    -0.91052213707850282454, -0.88931544599511413995, -0.86599939815409276989, -0.84062929625258031585,
    -0.81326531512279753855, -0.78397235894334138528, -0.75281990726053193974, -0.71988185017161077095,
    -0.68523631305423327031, -0.64896547125465731121, -0.61115535517239327756, -0.57189564620263400041,
    -0.53127946401989456504, -0.48940314570705295560, -0.44636601725346408687, -0.40227015796399162584,
    -0.35722015833766812554, -0.31132287199021096979, -0.26468716220876742362, -0.21742364374000708316,
    -0.16964442042399280330, -0.12146281929612055828, -0.07299312178779904237, -0.02435029266342442905,
    0.02435029266342442905, 0.07299312178779904237, 0.12146281929612055828, 0.16964442042399280330,
    0.21742364374000708316, 0.26468716220876742362, 0.31132287199021096979, 0.35722015833766812554,
    0.40227015796399162584, 0.44636601725346408687, 0.48940314570705295560, 0.53127946401989456504,
    0.57189564620263400041, 0.61115535517239327756, 0.64896547125465731121, 0.68523631305423327031,
    0.71988185017161077095, 0.75281990726053193974, 0.78397235894334138528, 0.81326531512279753855,
    0.84062929625258031585, 0.86599939815409276989, 0.88931544599511413995, 0.91052213707850282454,
    0.92956917213193956950, 0.94641137485840276522, 0.96100879965205376898, 0.97332682778991097550,
    0.98333625388462597705, 0.99101337147674428696, 0.99634011677195521983, 0.99930504173577217042};
static constexpr double GL64_W[64] = {
    0.00178328072169421517, 0.00414703326056292329, 0.00650445796897965427, 0.00884675982636439102,
    0.01116813946013146645, 0.01346304789671823147, 0.01572603047602508242, 0.01795171577569730156,
    0.02013482315353009450, 0.02227017380838300711, 0.02435270256871085309, 0.02637746971505462723,
    0.02833967261425970191, 0.03023465707240249531, 0.03205792835485145320, 0.03380516183714178668,
    0.03547221325688232341, 0.03705512854024015090, 0.03855015317861559127, 0.03995374113272034955,
    0.04126256324262348590, 0.04247351512365359766, 0.04358372452932346430, 0.04459055816375654541,
    0.04549162792741811429, 0.04628479658131437469, 0.04696818281620999957, 0.04754016571483030140,
    0.04799938859645831724, 0.04834476223480295431, 0.04857546744150345597, 0.04869095700913975144,
    0.04869095700913975144, 0.04857546744150345597, 0.04834476223480295431, 0.04799938859645831724,
    0.04754016571483030140, 0.04696818281620999957, 0.04628479658131437469, 0.04549162792741811429,
    0.04459055816375654541, 0.04358372452932346430, 0.04247351512365359766, 0.04126256324262348590,
    0.03995374113272034955, 0.03855015317861559127, 0.03705512854024015090, 0.03547221325688232341,
    0.03380516183714178668, 0.03205792835485145320, 0.03023465707240249531, 0.02833967261425970191,
    0.02637746971505462723, 0.02435270256871085309, 0.02227017380838300711, 0.02013482315353009450,
    0.01795171577569730156, 0.01572603047602508242, 0.01346304789671823147, 0.01116813946013146645,
    0.00884675982636439102, 0.00650445796897965427, 0.00414703326056292329, 0.00178328072169421517};

inline double smooth_step_erf(double E, double Ethr, double sigmaE) {
    if (sigmaE <= 0.0)
        return (E >= Ethr) ? 1.0 : 0.0;
    const double t = (E - Ethr) / (std::sqrt(2.0) * sigmaE);

    return std::max(0.0, std::min(1.0, 0.5 * (1.0 + std::erf(t))));
}

inline double Akin_smooth(double p,
                          double pthr_p, double pthr_pi,
                          double P, double alpha,
                          double sigmaE_p, double sigmaE_pi) {
    const double beta = beta_from_p(p);
    const double gamma = gamma_from_p(p);
    if (beta <= 0.0)
        return 0.0;

    const double lam = kallen(mL * mL, mp * mp, mpi * mpi);
    if (lam <= 0.0)
        return 0.0;

    const double pst = 0.5 * std::sqrt(lam) / mL;
    const double Epst = std::sqrt(mp * mp + pst * pst);
    const double Epist = std::sqrt(mpi * mpi + pst * pst);

    const double Ethr_p = std::sqrt(mp * mp + pthr_p * pthr_p);
    const double Ethr_pi = std::sqrt(mpi * mpi + pthr_pi * pthr_pi);

    const double A_p = gamma * Epst;
    const double B_p = gamma * beta * pst;
    const double A_pi = gamma * Epist;
    const double B_pi = -gamma * beta * pst;

    auto eps = [&](double x) -> double {
        const double Ep = A_p + B_p * x;
        const double Epi = A_pi + B_pi * x;
        const double ep = smooth_step_erf(Ep, Ethr_p, sigmaE_p);
        const double epi = smooth_step_erf(Epi, Ethr_pi, sigmaE_pi);
        return ep * epi;
    };

    double I0 = 0.0, I1 = 0.0;
    for (int i = 0; i < 64; ++i) {
        const double x = GL64_X[i];
        const double w = GL64_W[i];
        const double e = eps(x);
        I0 += w * e;
        I1 += w * (x * e);
    }

    const double A = 0.5 * I0 + 0.5 * alpha * P * I1;
    return A;
}

inline double Akin_with_cfg(double p, const AcceptCfg& C) {
    if (C.sigmaE_p > 0.0 || C.sigmaE_pi > 0.0) {
        return Akin_smooth(p, C.pthr_p, C.pthr_pi, C.PLambda, C.alpha,
                           C.sigmaE_p, C.sigmaE_pi);
    }
    return Akin(p, C.pthr_p, C.pthr_pi, C.PLambda, C.alpha);
}

inline double APS_from_p(double p) {
    const double beta_gamma = p / mL;
    const double Alen = Alength(beta_gamma, gCfg.Lmin_cm, gCfg.Lmax_cm);
    const double A2 = Akin_with_cfg(p, gCfg);
    const double aps = BpPi * Alen * A2;
    return aps;
}
inline double APS_from_bg(double beta_gamma) {
    const double p = beta_gamma * mL;
    return APS_from_p(p);
}

}

void draw_APS_vs_bg(const std::string& outdir) {
    const int N = 1200;
    std::vector<double> x(N), y_nom(N), y_loL(N), y_hiL(N), y_noThr(N), y_polp(N), y_polm(N);

    ana::AcceptCfg cfg_nom = ana::gCfg;
    ana::AcceptCfg cfg_loL = ana::gCfg;
    cfg_loL.Lmin_cm = 0.5;
    ana::AcceptCfg cfg_hiL = ana::gCfg;
    cfg_hiL.Lmin_cm = 2.5;
    ana::AcceptCfg cfg_thr = ana::gCfg;
    cfg_thr.pthr_p = 0.0;
    cfg_thr.pthr_pi = 0.0;
    cfg_thr.sigmaE_p = 0.0;
    cfg_thr.sigmaE_pi = 0.0;
    ana::AcceptCfg cfg_pP = ana::gCfg;
    cfg_pP.PLambda = +0.4;
    ana::AcceptCfg cfg_mP = ana::gCfg;
    cfg_mP.PLambda = -0.4;

    auto eval_APS = [](double bg, const ana::AcceptCfg& C) -> double {
        const double p = bg * ana::mL;
        const double Alen = ana::Alength(bg, C.Lmin_cm, C.Lmax_cm);
        const double A2 = ana::Akin_with_cfg(p, C);
        return ana::BpPi * Alen * A2;
    };

    const double bg_min = 0.0, bg_max = 10.0;
    for (int i = 0; i < N; ++i) {
        const double bg = bg_min + (bg_max - bg_min) * i / (N - 1);
        x[i] = bg;
        y_nom[i] = eval_APS(bg, cfg_nom);
        y_loL[i] = eval_APS(bg, cfg_loL);
        y_hiL[i] = eval_APS(bg, cfg_hiL);
        y_noThr[i] = eval_APS(bg, cfg_thr);
        y_polp[i] = eval_APS(bg, cfg_pP);
        y_polm[i] = eval_APS(bg, cfg_mP);
    }

    TCanvas c("c_APS", "APS vs beta*gamma", 900, 700);
    gStyle->SetOptStat(0);

    auto* pad_main = new TPad("pad_main", "pad_main", 0., 0.00, 1., 0.82);
    auto* pad_legend = new TPad("pad_legend", "pad_legend", 0., 0.82, 1., 1.00);
    pad_main->SetTopMargin(0.02);
    pad_main->SetBottomMargin(0.12);
    pad_main->SetLeftMargin(0.15);
    pad_main->SetRightMargin(0.05);
    pad_legend->SetTopMargin(0.05);
    pad_legend->SetBottomMargin(0.05);
    pad_legend->SetLeftMargin(0.02);
    pad_legend->SetRightMargin(0.02);
    pad_main->Draw();
    pad_legend->Draw();

    pad_main->cd();
    auto g_nom = new TGraph(N, x.data(), y_nom.data());
    auto g_loL = new TGraph(N, x.data(), y_loL.data());
    auto g_hiL = new TGraph(N, x.data(), y_hiL.data());
    auto g_thr = new TGraph(N, x.data(), y_noThr.data());
    auto g_polp = new TGraph(N, x.data(), y_polp.data());
    auto g_polm = new TGraph(N, x.data(), y_polm.data());

    g_nom->SetLineColor(kBlue + 1);
    g_nom->SetLineWidth(3);
    g_loL->SetLineColor(kGreen + 2);
    g_loL->SetLineStyle(2);
    g_loL->SetLineWidth(2);
    g_hiL->SetLineColor(kGreen + 2);
    g_hiL->SetLineStyle(9);
    g_hiL->SetLineWidth(2);
    g_thr->SetLineColor(kRed + 1);
    g_thr->SetLineStyle(3);
    g_thr->SetLineWidth(2);
    g_polp->SetLineColor(kMagenta + 1);
    g_polp->SetLineStyle(7);
    g_polp->SetLineWidth(2);
    g_polm->SetLineColor(kMagenta + 1);
    g_polm->SetLineStyle(7);
    g_polm->SetLineWidth(2);

    g_nom->SetTitle(";#beta#gamma ( #it{#Lambda} boost );A_{PS}(#beta#gamma)");
    g_nom->GetXaxis()->SetLimits(bg_min, bg_max);
    g_nom->GetYaxis()->SetRangeUser(0.0, 0.8);

    g_nom->Draw("AL");
    g_thr->Draw("L SAME");
    g_loL->Draw("L SAME");
    g_hiL->Draw("L SAME");
    g_polp->Draw("L SAME");
    g_polm->Draw("L SAME");

    pad_legend->cd();
    TLegend leg(0.12, 0.05, 0.95, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetNColumns(2);
    leg.AddEntry(g_nom, Form("Nominal thresholds (#sigma_{E}^{p}=%.0f MeV, #sigma_{E}^{#pi}=%.0f MeV)", 1000. * ana::gCfg.sigmaE_p, 1000. * ana::gCfg.sigmaE_pi), "l");
    leg.AddEntry(g_thr, "No thresholds", "l");
    leg.AddEntry(g_loL, "L_{min}=0.5 cm", "l");
    leg.AddEntry(g_hiL, "L_{min}=2.5 cm", "l");
    leg.AddEntry(g_polp, "P_{#Lambda}=+0.4", "l");
    leg.AddEntry(g_polm, "P_{#Lambda}=-0.4", "l");
    leg.Draw();

    c.SaveAs((outdir + "/APS_vs_bg.png").c_str());
    c.SaveAs((outdir + "/APS_vs_bg.pdf").c_str());
}

void draw_Akin_vs_bg(const std::string& outdir) {

    const int N = 1200;
    const double bg_min = 0.0, bg_max = 10.0;
    std::vector<double> x(N), y_nom(N), y_noThr(N), y_polp(N), y_polm(N);

    ana::AcceptCfg c_nom = ana::gCfg;
    c_nom.PLambda = 0.0;
    ana::AcceptCfg c_thr = ana::gCfg;
    c_thr.pthr_p = 0.0;
    c_thr.pthr_pi = 0.0;
    c_thr.sigmaE_p = 0.0;
    c_thr.sigmaE_pi = 0.0;
    ana::AcceptCfg c_pP = ana::gCfg;
    c_pP.PLambda = +0.4;
    ana::AcceptCfg c_mP = ana::gCfg;
    c_mP.PLambda = -0.4;

    for (int i = 0; i < N; ++i) {
        const double bg = bg_min + (bg_max - bg_min) * i / (N - 1);
        const double p = bg * ana::mL;
        x[i] = bg;
        y_nom[i] = ana::Akin_with_cfg(p, c_nom);
        y_noThr[i] = ana::Akin_with_cfg(p, c_thr);
        y_polp[i] = ana::Akin_with_cfg(p, c_pP);
        y_polm[i] = ana::Akin_with_cfg(p, c_mP);
    }

    TCanvas c("c_Akin", "A_kin vs beta*gamma", 900, 700);
    gStyle->SetOptStat(0);

    auto* pad_main = new TPad("pad_main", "pad_main", 0., 0.00, 1., 0.82);
    auto* pad_legend = new TPad("pad_legend", "pad_legend", 0., 0.82, 1., 1.00);
    pad_main->SetTopMargin(0.02);
    pad_main->SetBottomMargin(0.12);
    pad_main->SetLeftMargin(0.15);
    pad_main->SetRightMargin(0.05);
    pad_legend->SetTopMargin(0.05);
    pad_legend->SetBottomMargin(0.05);
    pad_legend->SetLeftMargin(0.02);
    pad_legend->SetRightMargin(0.02);
    pad_main->Draw();
    pad_legend->Draw();

    pad_main->cd();
    auto g_nom = new TGraph(N, x.data(), y_nom.data());
    auto g_thr = new TGraph(N, x.data(), y_noThr.data());
    auto g_polp = new TGraph(N, x.data(), y_polp.data());
    auto g_polm = new TGraph(N, x.data(), y_polm.data());

    g_nom->SetLineColor(kBlue + 1);
    g_nom->SetLineWidth(3);
    g_thr->SetLineColor(kRed + 1);
    g_thr->SetLineStyle(3);
    g_thr->SetLineWidth(2);
    g_polp->SetLineColor(kMagenta + 1);
    g_polp->SetLineStyle(7);
    g_polp->SetLineWidth(2);
    g_polm->SetLineColor(kMagenta + 1);
    g_polm->SetLineStyle(7);
    g_polm->SetLineWidth(2);

    g_nom->SetTitle(";#beta#gamma ( #it{#Lambda} boost );A_{kin}");
    g_nom->GetXaxis()->SetLimits(bg_min, bg_max);
    g_nom->GetYaxis()->SetRangeUser(0.0, 1.0);

    g_nom->Draw("AL");
    g_thr->Draw("L SAME");
    g_polp->Draw("L SAME");
    g_polm->Draw("L SAME");

    pad_legend->cd();
    TLegend leg(0.12, 0.05, 0.95, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetNColumns(2);
    leg.AddEntry(g_nom, Form("Nominal thresholds (#sigma_{E}^{p}=%.0f MeV, #sigma_{E}^{#pi}=%.0f MeV)", 1000. * ana::gCfg.sigmaE_p, 1000. * ana::gCfg.sigmaE_pi), "l");
    leg.AddEntry(g_thr, "No thresholds", "l");
    leg.AddEntry(g_polp, "P_{#Lambda}=+0.4", "l");
    leg.AddEntry(g_polm, "P_{#Lambda}=-0.4", "l");
    leg.Draw();

    c.SaveAs((outdir + "/Akin_vs_bg.png").c_str());
    c.SaveAs((outdir + "/Akin_vs_bg.pdf").c_str());
}

void draw_APS_heatmap_Lmin_vs_bg(const std::string& outdir) {

    const int Nx = 400, Ny = 120;
    const double bg_min = 0.0, bg_max = 10.0;
    const double Lmin_min = 0.0, Lmin_max = 3.0;

    TH2D h("hAPS", ";#beta#gamma ( #it{#Lambda} boost );L_{min} [cm]",
           Nx, bg_min, bg_max, Ny, Lmin_min, Lmin_max);

    ana::AcceptCfg cfg = ana::gCfg;

    for (int ix = 1; ix <= Nx; ++ix) {
        const double bg = h.GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= Ny; ++iy) {
            cfg.Lmin_cm = h.GetYaxis()->GetBinCenter(iy);
            const double p = bg * ana::mL;
            const double Alen = ana::Alength(bg, cfg.Lmin_cm, cfg.Lmax_cm);
            const double A2 = ana::Akin_with_cfg(p, cfg);
            const double aps = ana::BpPi * Alen * A2;
            h.SetBinContent(ix, iy, aps);
        }
    }

    TCanvas c("c_heat", "A_{PS}(#beta#gamma, L_{min})", 900, 700);
    gStyle->SetOptStat(0);
    c.SetTopMargin(0.08);
    c.SetBottomMargin(0.12);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.18);

    h.SetContour(60);
    h.GetZaxis()->SetRangeUser(0.0, 0.7);
    h.Draw("COLZ");
    c.Update();

    if (auto* palette = dynamic_cast<TPaletteAxis*>(h.GetListOfFunctions()->FindObject("palette"))) {
        palette->SetX1NDC(0.86);
        palette->SetX2NDC(0.90);
        palette->SetY1NDC(0.15);
        palette->SetY2NDC(0.92);
    }

    auto* h_contours = static_cast<TH2D*>(h.Clone("hAPS_contours"));
    if (h_contours) {
        h_contours->SetDirectory(nullptr);
        h_contours->SetLineColor(kGray + 3);
        h_contours->SetLineWidth(1);
        h_contours->Draw("CONT3 SAME");
    }
    c.SaveAs((outdir + "/APS_heatmap_Lmin_vs_bg.png").c_str());
    c.SaveAs((outdir + "/APS_heatmap_Lmin_vs_bg.pdf").c_str());
}

void analytical_correction() {

    ana::gCfg = {
        .Lmin_cm = 1.5,
        .Lmax_cm = 200.,
        .pthr_p = 0.25,
        .pthr_pi = 0.10,
        .PLambda = 0.0,
        .alpha = ana::alphaLambda,

        .sigmaE_p = 0.010,
        .sigmaE_pi = 0.010};

    const std::string outdir = "plots/analytical";
    gSystem->mkdir(outdir.c_str(), kTRUE);
    rarexsec::plot::Plotter{}.set_global_style();
    gStyle->SetNumberContours(60);

    std::cout << "[analytical] Using: Lmin=" << ana::gCfg.Lmin_cm << " cm, "
              << "Lmax=" << ana::gCfg.Lmax_cm << " cm, "
              << "p_thr^p=" << ana::gCfg.pthr_p << " GeV, "
              << "p_thr^pi=" << ana::gCfg.pthr_pi << " GeV, "
              << "P_Lambda=" << ana::gCfg.PLambda << ", "
              << "sigmaE_p=" << ana::gCfg.sigmaE_p << " GeV, "
              << "sigmaE_pi=" << ana::gCfg.sigmaE_pi << " GeV" << std::endl;

    draw_APS_vs_bg(outdir);
    draw_Akin_vs_bg(outdir);
    draw_APS_heatmap_Lmin_vs_bg(outdir);
}