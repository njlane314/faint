#include "TSystem.h"
#include "rarexsec/Hub.hh"
#include "rarexsec/proc/Selection.hh"
#include <iostream>
#include <string>
#include <vector>

namespace study {
inline bool is_numu_cc(int ch) { switch (ch) { case 10: case 11: case 12: case 13: case 15: case 16: case 18: return true; default: return false; } }
}

void main() {
  ROOT::EnableImplicitMT();
  const char* cfg  = gSystem->Getenv("RAREXSEC_CONFIG");
  const char* bl   = gSystem->Getenv("RAREXSEC_BEAMLINE");
  const char* pers = gSystem->Getenv("RAREXSEC_PERIODS");
  if (!cfg || !bl || !pers) { std::cerr << "missing env\n"; return; }
  rarexsec::Hub hub(cfg);
  auto periods = rarexsec::selection::split_csv(pers);
  auto mc = hub.simulation_entries(bl, periods);
  auto r = rarexsec::selection::evaluate(mc, study::is_numu_cc, rarexsec::selection::Preset::InclusiveMuCC);
  std::cout << std::fixed;
  std::cout << "Efficiency: " << r.efficiency() << "\n";
  std::cout << "Purity    : " << r.purity() << "\n";
}
