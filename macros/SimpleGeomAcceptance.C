// SimpleGeomAcceptance.C
// Geometry-only, purely analytical acceptance with SIX-side treatment,
// size-faithful Z×X and Z×Y plots, and overlays for:
//   - Standard FV (SFV) from your definition,
//   - Uniform-acceptance (UA) box from analytics,
//   - Intersection (IX) = inner box of SFV ∩ UA.
//
// Run:
//   root -l -q SimpleGeomAcceptance.C++
// or
//   root -l; .L SimpleGeomAcceptance.C+; SimpleGeomAcceptance();

#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TBox.h>
#include <TLegend.h>
#include <TMath.h>
#include <TString.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// ------------------ Hard-coded detector geometry (active TPC), cm ------------------
// MicroBooNE (active): x ∈ [0, 256.35], y ∈ [-116.5, 116.5], z ∈ [0, 1036.8]
static const double kXmin = 0.0,    kXmax = 256.35;
static const double kYmin = -116.5, kYmax = 116.5;
static const double kZmin = 0.0,    kZmax = 1036.8;

// ------------------ Standard fiducial volume (SFV) from your image, cm ------------------
static const double kSFV_Xmin =   5.0,  kSFV_Xmax = 251.0;
static const double kSFV_Ymin = -110.0, kSFV_Ymax = 110.0;
static const double kSFV_Zmin =  20.0,  kSFV_Zmax = 986.0;

// ------------------ Binning (for maps) ------------------
static const int kNX = 26, kNY = 24, kNZ = 52;

// ------------------ Analytic knobs (six-face, but symmetric defaults) ------------------
// Acceptance per face uses a simple exponential approach to plateau.
// Face-specific "cushion" S (cm) and slope L (cm) — you can tune asymmetrically if needed.
static const double kS_xmin = 45.0, kS_xmax = 45.0;
static const double kS_ymin = 45.0, kS_ymax = 45.0;
static const double kS_zmin = 45.0, kS_zmax = 45.0;
// Slopes (≈ βγ·cτ_Λ). Keep them face-specific in case you want asymmetries.
static const double kL_xmin = 3.945, kL_xmax = 3.945;
static const double kL_ymin = 3.945, kL_ymax = 3.945;
static const double kL_zmin = 3.945, kL_zmax = 3.945;

// Target uniform acceptance "inside the UA box" (plateau fraction)
static const double kA_uniform = 0.95;

// ------------------ Plot look ------------------
static const bool  kMakePNGs      = true;
static const char* kOutRoot       = "SimpleGeomAcceptance.root";
static const int   kBaseCanvasH   = 900;
static const double kLeftMargin   = 0.12, kRightMargin = 0.18, kTopMargin = 0.06, kBotMargin = 0.12;

// Box styles
static const int kColUA  = kRed+1,   kWUA  = 4;
static const int kColSFV = kAzure+2, kWSFV = 3;
static const int kColIX  = kGreen+2, kWIX  = 5, kLSIX = 7; // dashed

// Optional iso-contour at A = kA_uniform
static const bool kDrawContour  = true;
static const int  kContourColor = kBlack;

// ------------------ helpers ------------------
static inline double posPart(double x){ return (x>0.0)? x : 0.0; }

// Single-face acceptance with exponential approach:
// A_face(d; S,L) = 0                    if (d <= S)
//                = 1 - exp(-(d - S)/L) otherwise
static inline double A_face(double d_to_face, double S, double L){
  const double t = d_to_face - S;
  if (t <= 0.0) return 0.0;
  return 1.0 - std::exp(-t / L);
}

// Two-sided axis acceptance (averaged over ± direction along the axis):
// A_axis = 0.5*A_face(min-side) + 0.5*A_face(max-side)
static inline double A_axis_two_sided(double d_min, double d_max,
                                      double S_min, double S_max,
                                      double L_min, double L_max){
  const double a1 = A_face(d_min, S_min, L_min);
  const double a2 = A_face(d_max, S_max, L_max);
  return 0.5*(a1 + a2);
}

// Full event acceptance (geometry-only), product over axes
static inline double A_event(double x, double y, double z){
  const double dx_min = x - kXmin, dx_max = kXmax - x;
  const double dy_min = y - kYmin, dy_max = kYmax - y;
  const double dz_min = z - kZmin, dz_max = kZmax - z;

  const double Ax = A_axis_two_sided(dx_min, dx_max, kS_xmin, kS_xmax, kL_xmin, kL_xmax);
  const double Ay = A_axis_two_sided(dy_min, dy_max, kS_ymin, kS_ymax, kL_ymin, kL_ymax);
  const double Az = A_axis_two_sided(dz_min, dz_max, kS_zmin, kS_zmax, kL_zmin, kL_zmax);
  return Ax * Ay * Az;
}

static TCanvas* MakeSizeFaithfulCanvas(const char* name,
                                       double xrange, double yrange,
                                       int baseH = kBaseCanvasH,
                                       double lm=kLeftMargin, double rm=kRightMargin,
                                       double tm=kTopMargin,  double bm=kBotMargin) {
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

// Solve per-face margin to guarantee per-axis plateau fraction Aaxis per face
// (conservative: if each face meets Aaxis, then two-sided average ≥ Aaxis).
static inline double MarginFromTarget(double Aaxis, double S, double L){
  // A_face = 1 - exp(-(M - S)/L) >= Aaxis  ->  M = S - L*ln(1 - Aaxis)
  return S - L * std::log(1.0 - Aaxis);
}

// ------------------ main ------------------
void SimpleGeomAcceptance() {
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(100);

  // Target per-axis face fraction (conservative construction)
  const double Aaxis = std::pow(kA_uniform, 1.0/3.0);

  // Compute SIX per-face margins
  const double M_xmin = MarginFromTarget(Aaxis, kS_xmin, kL_xmin);
  const double M_xmax = MarginFromTarget(Aaxis, kS_xmax, kL_xmax);
  const double M_ymin = MarginFromTarget(Aaxis, kS_ymin, kL_ymin);
  const double M_ymax = MarginFromTarget(Aaxis, kS_ymax, kL_ymax);
  const double M_zmin = MarginFromTarget(Aaxis, kS_zmin, kL_zmin);
  const double M_zmax = MarginFromTarget(Aaxis, kS_zmax, kL_zmax);

  // Uniform-acceptance (UA) box (six-sided, not forced symmetric)
  const double UA_X1 = kXmin + M_xmin, UA_X2 = kXmax - M_xmax;
  const double UA_Y1 = kYmin + M_ymin, UA_Y2 = kYmax - M_ymax;
  const double UA_Z1 = kZmin + M_zmin, UA_Z2 = kZmax - M_zmax;

  // Intersection box = inner(SFV ∩ UA)
  const double IX_X1 = std::max(UA_X1, kSFV_Xmin);
  const double IX_X2 = std::min(UA_X2, kSFV_Xmax);
  const double IX_Y1 = std::max(UA_Y1, kSFV_Ymin);
  const double IX_Y2 = std::min(UA_Y2, kSFV_Ymax);
  const double IX_Z1 = std::max(UA_Z1, kSFV_Zmin);
  const double IX_Z2 = std::min(UA_Z2, kSFV_Zmax);
  const bool   kHasIX = (IX_X1<IX_X2) && (IX_Y1<IX_Y2) && (IX_Z1<IX_Z2);

  // Histograms with Z on x-axis:
  //   hZX: Z (x-axis) vs X (y-axis), averaged over Y
  //   hZY: Z (x-axis) vs Y (y-axis), averaged over X
  TH2D hZX("Acc_ZX","Predicted (geometry-only) acceptance;Z [cm];X [cm]",
           kNZ, kZmin, kZmax, kNX, kXmin, kXmax);
  TH2D hZY("Acc_ZY","Predicted (geometry-only) acceptance;Z [cm];Y [cm]",
           kNZ, kZmin, kZmax, kNY, kYmin, kYmax);

  // Fill Z×X (avg over Y)
  for (int iz=1; iz<=kNZ; ++iz) {
    const double z = hZX.GetXaxis()->GetBinCenter(iz);
    for (int ix=1; ix<=kNX; ++ix) {
      const double x = hZX.GetYaxis()->GetBinCenter(ix);
      double sum = 0.0;
      for (int iy=1; iy<=kNY; ++iy) {
        const double y = hZY.GetYaxis()->GetBinCenter(iy);
        sum += A_event(x,y,z);
      }
      hZX.SetBinContent(iz, ix, sum / double(kNY));
    }
  }

  // Fill Z×Y (avg over X)
  for (int iz=1; iz<=kNZ; ++iz) {
    const double z = hZY.GetXaxis()->GetBinCenter(iz);
    for (int iy=1; iy<=kNY; ++iy) {
      const double y = hZY.GetYaxis()->GetBinCenter(iy);
      double sum = 0.0;
      for (int ix=1; ix<=kNX; ++ix) {
        const double x = hZX.GetYaxis()->GetBinCenter(ix);
        sum += A_event(x,y,z);
      }
      hZY.SetBinContent(iz, iy, sum / double(kNX));
    }
  }

  // Save histos
  TFile fout(kOutRoot, "RECREATE");
  hZX.Write(); hZY.Write();
  fout.Close();

  // ------------- Plot: Z×X (Z on x-axis; size-faithful) -------------
  {
    TCanvas* c = MakeSizeFaithfulCanvas("c_ZX", kZmax-kZmin, kXmax-kXmin);
    hZX.SetMinimum(0); hZX.SetMaximum(1);
    hZX.Draw("COLZ");
    if (kDrawContour) {
      TH2D* hC = (TH2D*)hZX.Clone("hC_ZX");
      Double_t lvl[1] = {kA_uniform};
      hC->SetContour(1, lvl);
      hC->SetLineColor(kContourColor); hC->SetLineWidth(2);
      hC->Draw("CONT3 SAME");
    }
    // Boxes in (Z, X)
    TBox* bUA  = DrawBox(UA_Z1,  UA_X1,  UA_Z2,  UA_X2,  kColUA,  kWUA);
    TBox* bSFV = DrawBox(kSFV_Zmin,kSFV_Xmin,kSFV_Zmax,kSFV_Xmax,kColSFV,kWSFV);
    TBox* bIX  = nullptr;
    if (kHasIX) bIX = DrawBox(IX_Z1, IX_X1, IX_Z2, IX_X2, kColIX, kWIX, kLSIX);

    auto* leg = new TLegend(0.15, 0.90, 0.78, 0.99);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    leg->AddEntry(bUA,  Form("Uniform A #geq %.0f%% box", 100*kA_uniform), "l");
    leg->AddEntry(bSFV, "Standard FV box", "l");
    if (bIX) leg->AddEntry(bIX, "Intersection (inner) box", "l");
    leg->Draw();

    if (kMakePNGs) c->SaveAs("geom_acc_ZX.png");
  }

  // ------------- Plot: Z×Y (Z on x-axis; size-faithful) -------------
  {
    TCanvas* c = MakeSizeFaithfulCanvas("c_ZY", kZmax-kZmin, kYmax-kYmin);
    hZY.SetMinimum(0); hZY.SetMaximum(1);
    hZY.Draw("COLZ");
    if (kDrawContour) {
      TH2D* hC = (TH2D*)hZY.Clone("hC_ZY");
      Double_t lvl[1] = {kA_uniform};
      hC->SetContour(1, lvl);
      hC->SetLineColor(kContourColor); hC->SetLineWidth(2);
      hC->Draw("CONT3 SAME");
    }
    // Boxes in (Z, Y)
    TBox* bUA  = DrawBox(UA_Z1,  UA_Y1,  UA_Z2,  UA_Y2,  kColUA,  kWUA);
    TBox* bSFV = DrawBox(kSFV_Zmin,kSFV_Ymin,kSFV_Zmax,kSFV_Ymax,kColSFV,kWSFV);
    TBox* bIX  = nullptr;
    if (kHasIX) bIX = DrawBox(IX_Z1, IX_Y1, IX_Z2, IX_Y2, kColIX, kWIX, kLSIX);

    auto* leg = new TLegend(0.15, 0.90, 0.78, 0.99);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    leg->AddEntry(bUA,  Form("Uniform A #geq %.0f%% box", 100*kA_uniform), "l");
    leg->AddEntry(bSFV, "Standard FV box", "l");
    if (bIX) leg->AddEntry(bIX, "Intersection (inner) box", "l");
    leg->Draw();

    if (kMakePNGs) c->SaveAs("geom_acc_ZY.png");
  }

  // Console summary
  auto span = [](double a,double b){ return std::max(0.0, b-a); };
  const double V_UA = span(UA_X1,UA_X2)*span(UA_Y1,UA_Y2)*span(UA_Z1,UA_Z2);
  const double V_SF = span(kSFV_Xmin,kSFV_Xmax)*span(kSFV_Ymin,kSFV_Ymax)*span(kSFV_Zmin,kSFV_Zmax);
  const double V_IX = kHasIX ? span(IX_X1,IX_X2)*span(IX_Y1,IX_Y2)*span(IX_Z1,IX_Z2) : 0.0;

  std::cout << "\n=== Six-face per-axis margins to reach A >= " << 100*kA_uniform << "% (conservative) ===\n";
  std::cout << "Mx_min="<<M_xmin<<"  Mx_max="<<M_xmax
            << "  My_min="<<M_ymin<<"  My_max="<<M_ymax
            << "  Mz_min="<<M_zmin<<"  Mz_max="<<M_zmax<<"  [cm]\n";
  std::cout << "\nUA box (X:["<<UA_X1<<","<<UA_X2
            <<"], Y:["<<UA_Y1<<","<<UA_Y2
            <<"], Z:["<<UA_Z1<<","<<UA_Z2<<"])  Vol="<<V_UA<<" cm^3\n";
  std::cout << "SFV    (X:["<<kSFV_Xmin<<","<<kSFV_Xmax
            <<"], Y:["<<kSFV_Ymin<<","<<kSFV_Ymax
            <<"], Z:["<<kSFV_Zmin<<","<<kSFV_Zmax<<"])  Vol="<<V_SF<<" cm^3\n";
  if (kHasIX)
    std::cout << "IX     (X:["<<IX_X1<<","<<IX_X2
              <<"], Y:["<<IX_Y1<<","<<IX_Y2
              <<"], Z:["<<IX_Z1<<","<<IX_Z2<<"])  Vol="<<V_IX<<" cm^3\n";
  else
    std::cout << "IX     : (none — UA and SFV do not overlap along at least one axis)\n";

  std::cout << "Wrote " << kOutRoot << " and PNGs (if enabled).\n";
}

#ifndef __CINT__
int main(){ SimpleGeomAcceptance(); return 0; }
#endif
