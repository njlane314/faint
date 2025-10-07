// DeadWireImpact.C
// Analytical dead-wire impact with hard-coded configuration.
// Run:  root -l -q DeadWireImpact.C++
// Or:   root -l; .L DeadWireImpact.C+; DeadWireImpact();

#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

// ------------------ Hard-coded configuration ------------------
// Choose one mode:
static const bool kUseDeadWireLists = false;  // false = global fractions; true = read lists

// Global fractions (used if kUseDeadWireLists == false)
static const double kfU = 0.05;  // 5% dead on U
static const double kfV = 0.03;  // 3% dead on V
static const double kfY = 0.02;  // 2% dead on Y

// Event topology
static const int kNtracks = 3;   // mu, p, pi all need ≥2 planes

// Geometry for YZ maps (independent of X)
static const double kYmin = -116.5, kYmax = 116.5;
static const double kZmin =    0.0, kZmax = 1036.8;
static const int    kNY   = 24,     kNZ   = 52;

// Plane wiring (approx MicroBooNE-like)
struct Plane {
  const char* name; double theta_deg; double pitch_cm; int nWires; double wire0_offset;
};
static const Plane kU = {"U", +60.0, 0.3, 2400, 0.0};
static const Plane kV = {"V", -60.0, 0.3, 2400, 0.0};
static const Plane kY = {"Y",   0.0, 0.3, 2400, 0.0};

// Dead-wire lists (used if kUseDeadWireLists == true)
static const char* kDeadU = "dead_U.txt";
static const char* kDeadV = "dead_V.txt";
static const char* kDeadY = "dead_Y.txt";
// Local window (± wires) around vertex wire to estimate local live fraction
static const int kHalfWindow = 5;

// Output
static const char* kOutRoot  = "DeadWireImpact.root";
static const bool  kMakePNGs = true;
// --------------------------------------------------------------

struct DeadMask {
  int n=0; std::vector<char> live; std::vector<int> pref;
  void init(int nW){ n=nW; live.assign(n,1); }
  void markDead(int a,int b){
    if (n<=0) return;
    a=std::max(0,std::min(n-1,a)); b=std::max(0,std::min(n-1,b));
    if (a>b) std::swap(a,b);
    for (int i=a;i<=b;++i) live[i]=0;
  }
  void finalize(){
    pref.assign(n+1,0);
    for (int i=0;i<n;++i) pref[i+1]=pref[i]+(live[i]?1:0);
  }
  double fracLive(int L,int R) const{
    if (n<=0) return 1.0;
    L=std::max(0,std::min(n-1,L)); R=std::max(0,std::min(n-1,R));
    if (L>R) std::swap(L,R);
    const int cnt = R-L+1;
    const int alive = pref[R+1]-pref[L];
    return (cnt>0)? double(alive)/double(cnt) : 1.0;
  }
  double deadFrac() const { return (n>0)? 1.0 - double(pref.back())/double(n) : 0.0; }
};

static inline double deg2rad(double d){ return d * TMath::Pi() / 180.0; }
static inline int wireIndexYZ(const Plane& P, double y, double z){
  const double th = deg2rad(P.theta_deg);
  const double s  =  z*std::cos(th) - y*std::sin(th);
  int idx = int(std::floor(s/P.pitch_cm + P.wire0_offset + 0.5));
  if (idx<0) idx=0; if (idx>=P.nWires) idx=P.nWires-1;
  return idx;
}
static inline void loadRanges(const char* fn, DeadMask& m){
  if (!fn || !*fn) return;
  std::ifstream fin(fn);
  if (!fin) { std::cerr<<"[warn] cannot open "<<fn<<"; assuming none.\n"; return; }
  std::string line; int a,b;
  while (std::getline(fin,line)) {
    if (line.empty() || line[0]=='#') continue;
    std::istringstream iss(line); if (iss>>a>>b) m.markDead(a,b);
  }
}
static inline double R2_from_e(double eU,double eV,double eY){
  // P(track seen on >=2 planes) assuming independence across planes
  return eU*eV + eU*eY + eV*eY - 2.0*eU*eV*eY;
}

void DeadWireImpact() {
  gStyle->SetOptStat(0);

  // Output histos (YZ maps)
  TH2D hEU("eU_YZ","Plane U live fraction;Y [cm];Z [cm]", kNY, kYmin, kYmax, kNZ, kZmin, kZmax);
  TH2D hEV("eV_YZ","Plane V live fraction;Y [cm];Z [cm]", kNY, kYmin, kYmax, kNZ, kZmin, kZmax);
  TH2D hEY("eY_YZ","Plane Y live fraction;Y [cm];Z [cm]", kNY, kYmin, kYmax, kNZ, kZmin, kZmax);
  TH2D hR2("R2_YZ","Per-track P(#planes #geq 2);Y [cm];Z [cm]",
           kNY, kYmin, kYmax, kNZ, kZmin, kZmax);
  char ttl[128]; std::snprintf(ttl, sizeof(ttl), "Event readout factor (N=%d);Y [cm];Z [cm]", kNtracks);
  TH2D hEvt("Event_YZ", ttl, kNY, kYmin, kYmax, kNZ, kZmin, kZmax);

  double meanR2=0, meanEvt=0; int nb=0;

  if (!kUseDeadWireLists) {
    // -------- Global mode: constant maps from global fractions --------
    const double eU = 1.0 - kfU, eV = 1.0 - kfV, eY = 1.0 - kfY;
    const double R2 = R2_from_e(eU,eV,eY);
    const double Ev = std::pow(R2, std::max(1,kNtracks));

    for (int iy=1; iy<=kNY; ++iy)
      for (int iz=1; iz<=kNZ; ++iz) {
        hEU.SetBinContent(iy,iz,eU);
        hEV.SetBinContent(iy,iz,eV);
        hEY.SetBinContent(iy,iz,eY);
        hR2.SetBinContent(iy,iz,R2);
        hEvt.SetBinContent(iy,iz,Ev);
      }

    std::cout << "\n=== Dead-wire impact (GLOBAL) ===\n";
    std::cout << "fU="<<kfU<<", fV="<<kfV<<", fY="<<kfY
              << "  -> per-track R2="<<R2
              << " , event factor="<<Ev
              << "  (loss ~ "<<100.0*(1.0-Ev)<<"%)\n";
    meanR2=R2; meanEvt=Ev; nb=1;

  } else {
    // -------- List mode: read ranges, build local live-fraction maps --------
    DeadMask mU, mV, mY; mU.init(kU.nWires); mV.init(kV.nWires); mY.init(kY.nWires);
    loadRanges(kDeadU, mU); loadRanges(kDeadV, mV); loadRanges(kDeadY, mY);
    mU.finalize(); mV.finalize(); mY.finalize();

    for (int iy=1; iy<=kNY; ++iy) {
      const double y = hEU.GetXaxis()->GetBinCenter(iy);
      for (int iz=1; iz<=kNZ; ++iz) {
        const double z = hEU.GetYaxis()->GetBinCenter(iz);
        const int wU = wireIndexYZ(kU,y,z), wV = wireIndexYZ(kV,y,z), wY = wireIndexYZ(kY,y,z);
        const double eU = mU.fracLive(wU-kHalfWindow, wU+kHalfWindow);
        const double eV = mV.fracLive(wV-kHalfWindow, wV+kHalfWindow);
        const double eY = mY.fracLive(wY-kHalfWindow, wY+kHalfWindow);
        const double R2 = R2_from_e(eU,eV,eY);
        const double Ev = std::pow(R2, std::max(1,kNtracks));
        hEU.SetBinContent(iy,iz,eU);
        hEV.SetBinContent(iy,iz,eV);
        hEY.SetBinContent(iy,iz,eY);
        hR2.SetBinContent(iy,iz,R2);
        hEvt.SetBinContent(iy,iz,Ev);
        meanR2 += R2; meanEvt += Ev; nb++;
      }
    }

    meanR2 = (nb>0)? meanR2/nb : 1.0;
    meanEvt = (nb>0)? meanEvt/nb : 1.0;
    std::cout << "\n=== Dead-wire impact (LISTS) ===\n";
    std::cout << "Global dead fractions from masks: fU="<<mU.deadFrac()
              << ", fV="<<mV.deadFrac() << ", fY="<<mY.deadFrac() << "\n";
    std::cout << "Mean per-track R2 over YZ: " << meanR2
              << " ; mean event factor (N="<<kNtracks<<"): " << meanEvt
              << "  (loss ~ " << 100.0*(1.0-meanEvt) << "%)\n";
  }

  // Save + quick looks
  TFile fout(kOutRoot,"RECREATE");
  hEU.Write(); hEV.Write(); hEY.Write(); hR2.Write(); hEvt.Write();
  fout.Close();

  if (kMakePNGs) {
    TCanvas c("c","",1000,750);
    hEU.SetMinimum(0); hEU.SetMaximum(1); hEU.Draw("COLZ"); c.SaveAs("deadwire_eU_YZ.png");
    hEV.SetMinimum(0); hEV.SetMaximum(1); hEV.Draw("COLZ"); c.SaveAs("deadwire_eV_YZ.png");
    hEY.SetMinimum(0); hEY.SetMaximum(1); hEY.Draw("COLZ"); c.SaveAs("deadwire_eY_YZ.png");
    hR2.SetMinimum(0); hR2.SetMaximum(1); hR2.Draw("COLZ"); c.SaveAs("deadwire_R2_YZ.png");
    hEvt.SetMinimum(0); hEvt.SetMaximum(1); hEvt.Draw("COLZ"); c.SaveAs("deadwire_Event_YZ.png");
  }

  std::cout << "Wrote " << kOutRoot << " and PNGs (if enabled).\n";
}

#ifndef __CINT__
int main(){ DeadWireImpact(); return 0; }
#endif
