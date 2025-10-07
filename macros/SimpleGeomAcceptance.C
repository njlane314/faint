// SimpleGeomAcceptance.C
// Purely analytical, geometry-only acceptance with hard-coded configuration.
// Z is plotted on the x-axis in both views (ZxX and ZxY).
// Overlays: (1) uniform-acceptance box, (2) standard FV box, (3) their intersection box.
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
#include <TMath.h>
#include <TString.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// ------------------ Hard-coded configuration ------------------
// Active volume (MicroBooNE-like), cm
static const double kXmin = 0.0,    kXmax = 256.35;
static const double kYmin = -116.5, kYmax = 116.5;
static const double kZmin = 0.0,    kZmax = 1036.8;

// Standard FV (from your image), cm
static const double kSFV_Xmin =   5.0, kSFV_Xmax = 251.0;
static const double kSFV_Ymin = -110.0, kSFV_Ymax = 110.0;
static const double kSFV_Zmin =  20.0, kSFV_Zmax = 986.0;

// Binning for 2D maps
static const int kNX = 26, kNY = 24, kNZ = 52;

// Analytic acceptance knobs
// S_cushion: daughter path budget (e.g., max contained range) [cm]
static const double kS_cushion = 45.0;
// Exponential approach scales per axis Lx,Ly,Lz [cm] (≈ beta*gamma * c*tau_Λ).
// Here: beta*gamma ≈ 0.5, c*tau_Λ ≈ 7.89 cm -> 3.945 cm. Tune if desired.
static const double kLx = 3.945, kLy = 3.945, kLz = 3.945;

// Uniform acceptance target (inside the bound box we draw)
static const double kA_uniform = 0.95;   // plateau fraction target

// Plot appearance
static const bool  kMakePNGs      = true;
static const char* kOutRoot       = "SimpleGeomAcceptance.root";
static const int   kBaseCanvasH   = 900;      // base canvas height (px)
static const double kLeftMargin   = 0.12;     // pad margins
static const double kRightMargin  = 0.18;     // leave room for palette
static const double kTopMargin    = 0.06;
static const double kBotMargin    = 0.12;

// Box styles
static const int kColUniform = kRed+1;      // uniform-acceptance box
static const int kColSFV     = kAzure+2;    // standard FV box
static const int kColIntersect= kGreen+2;   // intersection box
static const int kWUniform   = 4;
static const int kWSFV       = 3;
static const int kWIntersect = 5;
static const int kLSIntersect= 7;           // dashed

// Optional: draw iso-contour at A = kA_uniform
static const bool kDrawContour  = true;
static const int  kContourColor = kBlack;
// --------------------------------------------------------------

// ------- helpers -------
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
static TCanvas* MakeSizeFaithfulCanvas(const char* name,
                                       double xrange, double yrange,
                                       int baseH = kBaseCanvasH,
                                       double lm=kLeftMargin, double rm=kRightMargin,
                                       double tm=kTopMargin,  double bm=kBotMargin) {
  // Ensure the inner (frame) area has aspect ratio = xrange/yrange
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
static TBox* DrawBox(double x1,double y1,double x2,double y2,int color,int width,int lstyle=1) {
  auto* b = new TBox(x1,y1,x2,y2);
  b->SetFillStyle(0);
  b->SetLineColor(color);
  b->SetLineWidth(width);
  b->SetLineStyle(lstyle);
  b->Draw("same");
  return b;
}

// ------- main -------
void SimpleGeomAcceptance() {
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(100);

  // --- Uniform-acceptance margins (symmetric inset per axis) ---
  const double Aaxis = std::pow(kA_uniform, 1.0/3.0);
  auto M_of = [&](double L){ return kS_cushion - L * std::log(1.0 - Aaxis); };
  const double Mx = M_of(kLx);
  const double My = M_of(kLy);
  const double Mz = M_of(kLz);

  // Uniform-acceptance box
  const double UA_X1 = kXmin + Mx, UA_X2 = kXmax - Mx;
  const double UA_Y1 = kYmin + My, UA_Y2 = kYmax - My;
  const double UA_Z1 = kZmin + Mz, UA_Z2 = kZmax - Mz;

  // Intersection box: inner box of (Uniform-acceptance) ∩ (Standard FV)
  const double IX_X1 = std::max(UA_X1, kSFV_Xmin);
  const double IX_X2 = std::min(UA_X2, kSFV_Xmax);
  const double IX_Y1 = std::max(UA_Y1, kSFV_Ymin);
  const double IX_Y2 = std::min(UA_Y2, kSFV_Ymax);
  const double IX_Z1 = std::max(UA_Z1, kSFV_Zmin);
  const double IX_Z2 = std::min(UA_Z2, kSFV_Zmax);
  const bool   kHasIntersection = (IX_X1<IX_X2) && (IX_Y1<IX_Y2) && (IX_Z1<IX_Z2);

  // Histograms where Z is on x-axis
  // hZX: Z (x-axis) vs X (y-axis) — averaged over Y
  TH2D hZX("Acc_ZX","Predicted (geometry-only) acceptance;Z [cm];X [cm]",
           kNZ, kZmin, kZmax, kNX, kXmin, kXmax);
  // hZY: Z (x-axis) vs Y (y-axis) — averaged over X
  TH2D hZY("Acc_ZY","Predicted (geometry-only) acceptance;Z [cm];Y [cm]",
           kNZ, kZmin, kZmax, kNY, kYmin, kYmax);

  // --- Fill Z×X (average over Y) ---
  for (int iz=1; iz<=kNZ; ++iz) {
    const double z = hZX.GetXaxis()->GetBinCenter(iz);
    for (int ix=1; ix<=kNX; ++ix) {
      const double x = hZX.GetYaxis()->GetBinCenter(ix);
      double sum=0;
      for (int iy=1; iy<=kNY; ++iy) {
        const double y = hZY.GetYaxis()->GetBinCenter(iy);
        const double dx = d_to_nearest(x, kXmin, kXmax);
        const double dy = d_to_nearest(y, kYmin, kYmax);
        const double dz = d_to_nearest(z, kZmin, kZmax);
        const double ax = axis_accept(dx, kS_cushion, kLx);
        const double ay = axis_accept(dy, kS_cushion, kLy);
        const double az = axis_accept(dz, kS_cushion, kLz);
        sum += ax*ay*az;
      }
      hZX.SetBinContent(iz, ix, sum / double(kNY));
    }
  }

  // --- Fill Z×Y (average over X) ---
  for (int iz=1; iz<=kNZ; ++iz) {
    const double z = hZY.GetXaxis()->GetBinCenter(iz);
    for (int iy=1; iy<=kNY; ++iy) {
      const double y = hZY.GetYaxis()->GetBinCenter(iy);
      double sum=0;
      for (int ix=1; ix<=kNX; ++ix) {
        const double x = hZX.GetYaxis()->GetBinCenter(ix);
        const double dx = d_to_nearest(x, kXmin, kXmax);
        const double dy = d_to_nearest(y, kYmin, kYmax);
        const double dz = d_to_nearest(z, kZmin, kZmax);
        const double ax = axis_accept(dx, kS_cushion, kLx);
        const double ay = axis_accept(dy, kS_cushion, kLy);
        const double az = axis_accept(dz, kS_cushion, kLz);
        sum += ax*ay*az;
      }
      hZY.SetBinContent(iz, iy, sum / double(kNX));
    }
  }

  // Save histograms
  TFile fout(kOutRoot, "RECREATE");
  hZX.Write(); hZY.Write();
  fout.Close();

  // -------------------- Plot: Z×X (size-faithful; Z on x-axis) --------------------
  {
    TCanvas* c = MakeSizeFaithfulCanvas("cZX", kZmax-kZmin, kXmax-kXmin);
    hZX.SetMinimum(0); hZX.SetMaximum(1);
    hZX.Draw("COLZ");
    if (kDrawContour) {
      TH2D* hC = (TH2D*)hZX.Clone("hC_ZX");
      Double_t lvl[1] = {kA_uniform};
      hC->SetContour(1, lvl);
      hC->SetLineColor(kContourColor); hC->SetLineWidth(2);
      hC->Draw("CONT3 SAME");
    }
    // Boxes on this view are in (Z, X)
    TBox* bUA = DrawBox(UA_Z1, UA_X1, UA_Z2, UA_X2, kColUniform, kWUniform);
    TBox* bSF = DrawBox(kSFV_Zmin, kSFV_Xmin, kSFV_Zmax, kSFV_Xmax, kColSFV, kWSFV);
    TBox* bIX = nullptr;
    if (kHasIntersection)
      bIX = DrawBox(IX_Z1, IX_X1, IX_Z2, IX_X2, kColIntersect, kWIntersect, kLSIntersect);

    // Legend
    auto* leg = new TLegend(0.15, 0.90, 0.70, 0.99);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);
    leg->AddEntry(bUA, Form("Uniform A #geq %.0f%% box", 100*kA_uniform), "l");
    leg->AddEntry(bSF, "Standard FV box", "l");
    if (bIX) leg->AddEntry(bIX, "Intersection (inner) box", "l");
    leg->Draw();

    if (kMakePNGs) c->SaveAs("geom_acc_ZX.png");
  }

  // -------------------- Plot: Z×Y (size-faithful; Z on x-axis) --------------------
  {
    TCanvas* c = MakeSizeFaithfulCanvas("cZY", kZmax-kZmin, kYmax-kYmin);
    hZY.SetMinimum(0); hZY.SetMaximum(1);
    hZY.Draw("COLZ");
    if (kDrawContour) {
      TH2D* hC = (TH2D*)hZY.Clone("hC_ZY");
      Double_t lvl[1] = {kA_uniform};
      hC->SetContour(1, lvl);
      hC->SetLineColor(kContourColor); hC->SetLineWidth(2);
      hC->Draw("CONT3 SAME");
    }
    // Boxes on this view are in (Z, Y)
    TBox* bUA = DrawBox(UA_Z1, UA_Y1, UA_Z2, UA_Y2, kColUniform, kWUniform);
    TBox* bSF = DrawBox(kSFV_Zmin, kSFV_Ymin, kSFV_Zmax, kSFV_Ymax, kColSFV, kWSFV);
    TBox* bIX = nullptr;
    if (kHasIntersection)
      bIX = DrawBox(IX_Z1, IX_Y1, IX_Z2, IX_Y2, kColIntersect, kWIntersect, kLSIntersect);

    // Legend
    auto* leg = new TLegend(0.15, 0.90, 0.70, 0.99);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);
    leg->AddEntry(bUA, Form("Uniform A #geq %.0f%% box", 100*kA_uniform), "l");
    leg->AddEntry(bSF, "Standard FV box", "l");
    if (bIX) leg->AddEntry(bIX, "Intersection (inner) box", "l");
    leg->Draw();

    if (kMakePNGs) c->SaveAs("geom_acc_ZY.png");
  }

  // Console summary
  auto span = [](double a,double b){ return std::max(0.0, b-a); };
  const double UA_Vol = span(UA_X1,UA_X2)*span(UA_Y1,UA_Y2)*span(UA_Z1,UA_Z2);
  const double SF_Vol = span(kSFV_Xmin,kSFV_Xmax)*span(kSFV_Ymin,kSFV_Ymax)*span(kSFV_Zmin,kSFV_Zmax);
  const double IX_Vol = kHasIntersection ? span(IX_X1,IX_X2)*span(IX_Y1,IX_Y2)*span(IX_Z1,IX_Z2) : 0.0;

  std::cout << "\n=== Boxes (cm) ===\n";
  std::cout << "Uniform-acceptance (A >= " << 100*kA_uniform << "%):\n"
            << "  X: ["<<UA_X1<<", "<<UA_X2<<"]  Y: ["<<UA_Y1<<", "<<UA_Y2<<"]  Z: ["<<UA_Z1<<", "<<UA_Z2<<"]\n";
  std::cout << "Standard FV:\n"
            << "  X: ["<<kSFV_Xmin<<", "<<kSFV_Xmax<<"]  Y: ["<<kSFV_Ymin<<", "<<kSFV_Ymax<<"]  Z: ["<<kSFV_Zmin<<", "<<kSFV_Zmax<<"]\n";
  if (kHasIntersection) {
    std::cout << "Intersection (inner) box:\n"
              << "  X: ["<<IX_X1<<", "<<IX_X2<<"]  Y: ["<<IX_Y1<<", "<<IX_Y2<<"]  Z: ["<<IX_Z1<<", "<<IX_Z2<<"]\n";
  } else {
    std::cout << "Intersection: (none — boxes do not overlap in at least one axis)\n";
  }
  std::cout << "\nVolumes (cm^3):\n"
            << "  Uniform-acceptance: " << UA_Vol << "\n"
            << "  Standard FV:        " << SF_Vol << "\n"
            << "  Intersection:       " << IX_Vol << "\n";
  std::cout << "Wrote " << kOutRoot << " and PNGs (if enabled).\n";
}

#ifndef __CINT__
int main(){ SimpleGeomAcceptance(); return 0; }
#endif
