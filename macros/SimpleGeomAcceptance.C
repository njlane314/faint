// SimpleGeomAcceptance.C
// Purely analytical, geometry-only acceptance with hard-coded configuration.
// Size-faithful plots + overlay of a "uniform acceptance" bounding box.
//
// Run (no args):
//   root -l -q SimpleGeomAcceptance.C++
// or
//   root -l; .L SimpleGeomAcceptance.C+; SimpleGeomAcceptance();

#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TBox.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>

// ------------------ Hard-coded configuration ------------------
// Active volume (MicroBooNE-like), cm
static const double kXmin = 0.0,    kXmax = 256.35;
static const double kYmin = -116.5, kYmax = 116.5;
static const double kZmin = 0.0,    kZmax = 1036.8;

// Binning for 2D maps
static const int kNX = 26, kNY = 24, kNZ = 52;

// Analytic acceptance knobs
// S_cushion: daughter path budget (e.g., max contained range) [cm]
static const double kS_cushion = 45.0;
// Exponential approach scales per axis Lx,Ly,Lz [cm] (≈ beta*gamma * c*tau_Λ).
// Here: beta*gamma ≈ 0.5, c*tau_Λ ≈ 7.89 cm -> 3.945 cm. Tune if desired.
static const double kLx = 3.945, kLy = 3.945, kLz = 3.945;

// Uniform acceptance target (inside the bound-box we draw)
static const double kA_uniform = 0.95;   // plateau fraction target

// Plot appearance
static const bool  kMakePNGs      = true;
static const char* kOutRoot       = "SimpleGeomAcceptance.root";
static const int   kBaseCanvasH   = 900;    // base canvas height (px); width auto-set for aspect
// Margins (fractions of pad). Right margin leaves space for the colorbar.
static const double kLeftMargin   = 0.12;
static const double kRightMargin  = 0.18;   // room for palette
static const double kTopMargin    = 0.06;
static const double kBotMargin    = 0.12;
static const int    kBoxColor     = kRed+1;
static const int    kContourColor = kBlack; // optional contour at A = kA_uniform
static const bool   kDrawContour  = true;
// --------------------------------------------------------------

static inline double axis_accept(double d, double S, double L) {
  const double t = d - S;
  if (t <= 0.0) return 0.0;
  return 1.0 - std::exp(-t / L);
}
static inline double d_to_nearest(double a, double amin, double amax) {
  const double left = a - amin;
  const double right = amax - a;
  return (left < right) ? left : right;
}
static inline double mean(const std::vector<double>& v) {
  if (v.empty()) return 0.0;
  double s=0; for (double x : v) s += x; return s / double(v.size());
}

static TCanvas* MakeSizeFaithfulCanvas(const char* name,
                                       double xrange, double yrange,
                                       int baseH = kBaseCanvasH,
                                       double lm=kLeftMargin, double rm=kRightMargin,
                                       double tm=kTopMargin,  double bm=kBotMargin) {
  // Ensure the inner (frame) area has aspect ratio = xrange/yrange
  // (W*(1-lm-rm)) / (H*(1-tm-bm)) = xrange/yrange  ->  W/H = (x/y)*((1-tm-bm)/(1-lm-rm))
  const double innerX = 1.0 - lm - rm;
  const double innerY = 1.0 - tm - bm;
  const double target = (xrange/yrange) * (innerY/innerX);
  const int H = baseH;
  const int W = std::max(600, int(std::round(target * H)));
  TCanvas* c = new TCanvas(name, name, W, H);
  c->SetLeftMargin(lm); c->SetRightMargin(rm);
  c->SetTopMargin(tm);  c->SetBottomMargin(bm);
  return c;
}
static void DrawBoundBox(double x1,double y1,double x2,double y2, int color=kBoxColor, int width=4) {
  auto* b = new TBox(x1,y1,x2,y2);
  b->SetFillStyle(0);
  b->SetLineColor(color);
  b->SetLineWidth(width);
  b->Draw("same");
}
static void LabelBoxDims(double x1,double y1,double x2,double y2, const char* coordLabel) {
  TLatex t;
  t.SetNDC(false);
  t.SetTextFont(42);
  t.SetTextSize(0.03);
  char buf[256];
  std::snprintf(buf,sizeof(buf), "%s box: [%.1f, %.1f] × [%.1f, %.1f] cm",
                coordLabel, x1, x2, y1, y2);
  // place near top-left corner inside box
  t.DrawLatex(x1 + 0.02*(x2-x1), y2 - 0.05*(y2-y1), buf);
}

void SimpleGeomAcceptance() {
  gStyle->SetOptStat(0);

  // Compute per-axis symmetric margins for the target acceptance
  // A_target = (Ax * Ay * Az) with Ax=1-exp(-(Mx-S)/Lx) at the bound
  const double Aaxis = std::pow(kA_uniform, 1.0/3.0);
  auto M_of = [&](double L){ return kS_cushion - L * std::log(1.0 - Aaxis); };
  const double Mx = M_of(kLx);
  const double My = M_of(kLy);
  const double Mz = M_of(kLz);

  // Bound-box corners for each view
  const double X1 = kXmin + Mx, X2 = kXmax - Mx;
  const double Y1 = kYmin + My, Y2 = kYmax - My;
  const double Z1 = kZmin + Mz, Z2 = kZmax - Mz;

  // Histograms (maps averaged over the suppressed axis)
  TH2D hXZ("Acc_XZ","Predicted (geometry-only) acceptance;X [cm];Z [cm]",
           kNX, kXmin, kXmax, kNZ, kZmin, kZmax);
  TH2D hXY("Acc_XY","Predicted (geometry-only) acceptance;X [cm];Y [cm]",
           kNX, kXmin, kXmax, kNY, kYmin, kYmax);
  TH2D hYZ("Acc_YZ","Predicted (geometry-only) acceptance;Y [cm];Z [cm]",
           kNY, kYmin, kYmax, kNZ, kZmin, kZmax);

  // Fill XZ (avg over Y)
  for (int ix=1; ix<=kNX; ++ix) {
    const double x = hXZ.GetXaxis()->GetBinCenter(ix);
    for (int iz=1; iz<=kNZ; ++iz) {
      const double z = hXZ.GetYaxis()->GetBinCenter(iz);
      std::vector<double> vals; vals.reserve(kNY);
      for (int iy=1; iy<=kNY; ++iy) {
        const double y = hXY.GetYaxis()->GetBinCenter(iy);
        const double dx = d_to_nearest(x, kXmin, kXmax);
        const double dy = d_to_nearest(y, kYmin, kYmax);
        const double dz = d_to_nearest(z, kZmin, kZmax);
        const double ax = axis_accept(dx, kS_cushion, kLx);
        const double ay = axis_accept(dy, kS_cushion, kLy);
        const double az = axis_accept(dz, kS_cushion, kLz);
        vals.push_back(ax*ay*az);
      }
      hXZ.SetBinContent(ix, iz, mean(vals));
    }
  }

  // Fill XY (avg over Z)
  for (int ix=1; ix<=kNX; ++ix) {
    const double x = hXY.GetXaxis()->GetBinCenter(ix);
    for (int iy=1; iy<=kNY; ++iy) {
      const double y = hXY.GetYaxis()->GetBinCenter(iy);
      std::vector<double> vals; vals.reserve(kNZ);
      for (int iz=1; iz<=kNZ; ++iz) {
        const double z = hXZ.GetYaxis()->GetBinCenter(iz);
        const double dx = d_to_nearest(x, kXmin, kXmax);
        const double dy = d_to_nearest(y, kYmin, kYmax);
        const double dz = d_to_nearest(z, kZmin, kZmax);
        const double ax = axis_accept(dx, kS_cushion, kLx);
        const double ay = axis_accept(dy, kS_cushion, kLy);
        const double az = axis_accept(dz, kS_cushion, kLz);
        vals.push_back(ax*ay*az);
      }
      hXY.SetBinContent(ix, iy, mean(vals));
    }
  }

  // Fill YZ (avg over X)
  for (int iy=1; iy<=kNY; ++iy) {
    const double y = hYZ.GetXaxis()->GetBinCenter(iy);
    for (int iz=1; iz<=kNZ; ++iz) {
      const double z = hYZ.GetYaxis()->GetBinCenter(iz);
      std::vector<double> vals; vals.reserve(kNX);
      for (int ix=1; ix<=kNX; ++ix) {
        const double x = hXY.GetXaxis()->GetBinCenter(ix);
        const double dx = d_to_nearest(x, kXmin, kXmax);
        const double dy = d_to_nearest(y, kYmin, kYmax);
        const double dz = d_to_nearest(z, kZmin, kZmax);
        const double ax = axis_accept(dx, kS_cushion, kLx);
        const double ay = axis_accept(dy, kS_cushion, kLy);
        const double az = axis_accept(dz, kS_cushion, kLz);
        vals.push_back(ax*ay*az);
      }
      hYZ.SetBinContent(iy, iz, mean(vals));
    }
  }

  // Save histograms
  TFile fout(kOutRoot, "RECREATE");
  hXZ.Write(); hXY.Write(); hYZ.Write();
  fout.Close();

  // -------------------- Plot: XZ (size-faithful) --------------------
  {
    TCanvas* c = MakeSizeFaithfulCanvas("cXZ", kXmax-kXmin, kZmax-kZmin);
    hXZ.SetMinimum(0); hXZ.SetMaximum(1);
    hXZ.Draw("COLZ");
    // Optional contour exactly at A = kA_uniform
    if (kDrawContour) {
      TH2D* hC = (TH2D*)hXZ.Clone("hC_XZ");
      double lvl = kA_uniform; hC->SetContour(1, &lvl);
      hC->SetLineColor(kContourColor); hC->SetLineWidth(2);
      hC->Draw("CONT3 SAME");
    }
    // Bound box overlay: X∈[X1,X2], Z∈[Z1,Z2]
    DrawBoundBox(X1, Z1, X2, Z2);
    LabelBoxDims(X1, Z1, X2, Z2, "XZ");
    // Legend
    auto* leg = new TLegend(0.15, 0.90, 0.60, 0.99);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);
    leg->AddEntry((TObject*)0, Form("Uniform A #geq %.0f%% box", 100*kA_uniform), "");
    leg->Draw();
    if (kMakePNGs) c->SaveAs("geom_acc_XZ.png");
  }

  // -------------------- Plot: XY (size-faithful) --------------------
  {
    TCanvas* c = MakeSizeFaithfulCanvas("cXY", kXmax-kXmin, kYmax-kYmin);
    hXY.SetMinimum(0); hXY.SetMaximum(1);
    hXY.Draw("COLZ");
    if (kDrawContour) {
      TH2D* hC = (TH2D*)hXY.Clone("hC_XY");
      double lvl = kA_uniform; hC->SetContour(1, &lvl);
      hC->SetLineColor(kContourColor); hC->SetLineWidth(2);
      hC->Draw("CONT3 SAME");
    }
    DrawBoundBox(X1, Y1, X2, Y2);
    LabelBoxDims(X1, Y1, X2, Y2, "XY");
    auto* leg = new TLegend(0.15, 0.90, 0.60, 0.99);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);
    leg->AddEntry((TObject*)0, Form("Uniform A #geq %.0f%% box", 100*kA_uniform), "");
    leg->Draw();
    if (kMakePNGs) c->SaveAs("geom_acc_XY.png");
  }

  // -------------------- Plot: YZ (size-faithful) --------------------
  {
    TCanvas* c = MakeSizeFaithfulCanvas("cYZ", kYmax-kYmin, kZmax-kZmin);
    hYZ.SetMinimum(0); hYZ.SetMaximum(1);
    hYZ.Draw("COLZ");
    if (kDrawContour) {
      TH2D* hC = (TH2D*)hYZ.Clone("hC_YZ");
      double lvl = kA_uniform; hC->SetContour(1, &lvl);
      hC->SetLineColor(kContourColor); hC->SetLineWidth(2);
      hC->Draw("CONT3 SAME");
    }
    DrawBoundBox(Y1, Z1, Y2, Z2);
    LabelBoxDims(Y1, Z1, Y2, Z2, "YZ");
    auto* leg = new TLegend(0.15, 0.90, 0.60, 0.99);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);
    leg->AddEntry((TObject*)0, Form("Uniform A #geq %.0f%% box", 100*kA_uniform), "");
    leg->Draw();
    if (kMakePNGs) c->SaveAs("geom_acc_YZ.png");
  }

  // Console summary
  std::cout << "\n=== Uniform-acceptance bound box (A >= " << 100*kA_uniform << "%) ===\n";
  std::cout << "Mx ~ " << Mx << " cm,  My ~ " << My << " cm,  Mz ~ " << Mz << " cm\n";
  std::cout << "X: [" << X1 << ", " << X2 << "] cm  (width ~ " << (X2-X1) << " cm)\n";
  std::cout << "Y: [" << Y1 << ", " << Y2 << "] cm  (height ~ " << (Y2-Y1) << " cm)\n";
  std::cout << "Z: [" << Z1 << ", " << Z2 << "] cm  (length ~ " << (Z2-Z1) << " cm)\n";
  std::cout << "Wrote " << kOutRoot << " and PNGs (if enabled).\n";
}

#ifndef __CINT__
int main(){ SimpleGeomAcceptance(); return 0; }
#endif
