// SimpleGeomAcceptance.C
// Purely analytical, geometry-only acceptance with hard-coded configuration.
// Run:  root -l -q SimpleGeomAcceptance.C++
// Or:   root -l; .L SimpleGeomAcceptance.C+; SimpleGeomAcceptance();

#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>

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
// Here: beta*gamma ≈ 0.5, c*tau_Λ ≈ 7.89 cm -> 3.945 cm. Set independently if desired.
static const double kLx = 3.945, kLy = 3.945, kLz = 3.945;

// Uniform acceptance target (everywhere inside a symmetric FV inset)
static const double kA_uniform = 0.95;

// Output
static const char* kOutRoot  = "SimpleGeomAcceptance.root";
static const bool  kMakePNGs = true;
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
  double s=0; for (double x : v) s += x;
  return s / double(v.size());
}

void SimpleGeomAcceptance() {
  gStyle->SetOptStat(0);

  // Histograms
  TH2D hXZ("Acc_XZ","Predicted (geometry-only) acceptance;X [cm];Z [cm]",
           kNX, kXmin, kXmax, kNZ, kZmin, kZmax);
  TH2D hXY("Acc_XY","Predicted (geometry-only) acceptance;X [cm];Y [cm]",
           kNX, kXmin, kXmax, kNY, kYmin, kYmax);
  TH2D hYZ("Acc_YZ","Predicted (geometry-only) acceptance;Y [cm];Z [cm]",
           kNY, kYmin, kYmax, kNZ, kZmin, kZmax);

  // --- XZ map (avg over Y)
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

  // --- XY map (avg over Z)
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

  // --- YZ map (avg over X)
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

  // Save + quick looks
  TFile fout(kOutRoot, "RECREATE");
  hXZ.Write(); hXY.Write(); hYZ.Write();
  fout.Close();

  if (kMakePNGs) {
    TCanvas c("c","",1000,750);
    hXZ.SetMinimum(0); hXZ.SetMaximum(1); hXZ.Draw("COLZ"); c.SaveAs("geom_acc_XZ.png");
    hXY.SetMinimum(0); hXY.SetMaximum(1); hXY.Draw("COLZ"); c.SaveAs("geom_acc_XY.png");
    hYZ.SetMinimum(0); hYZ.SetMaximum(1); hYZ.Draw("COLZ"); c.SaveAs("geom_acc_YZ.png");
  }

  // Uniform-FV margins (symmetric inset)
  if (kA_uniform>0 && kA_uniform<1) {
    const double Aaxis = std::pow(kA_uniform, 1.0/3.0);
    auto M = [&](double L){ return kS_cushion - L * std::log(1.0 - Aaxis); };
    std::cout << "\n=== Suggested symmetric margins for A >= " << kA_uniform*100 << "% ===\n";
    std::cout << "Mx ~ " << M(kLx) << " cm,  My ~ " << M(kLy) << " cm,  Mz ~ " << M(kLz) << " cm\n";
  }

  std::cout << "Wrote " << kOutRoot << " and PNGs (if enabled).\n";
}

#ifndef __CINT__
int main(){ SimpleGeomAcceptance(); return 0; }
#endif
