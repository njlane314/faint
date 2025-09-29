#include <iostream>
#include <memory>
#include <stdexcept>

#include "TSystem.h"

#include <rarexsec/MuonSelector.h>
#include <rarexsec/PreSelection.h>
#include <rarexsec/TruthClassifier.h>

void example_macro()
{
  if (gSystem->Load("librarexsec_root")) {
    throw std::runtime_error("Failed to load librexsec_root library");
  }

  auto preselection = std::make_unique<analysis::PreSelection>();
  auto muon = std::make_unique<analysis::MuonSelector>();
  auto truth = std::make_unique<analysis::TruthClassifier>();

  muon->chain_processor(std::move(truth));
  preselection->chain_processor(std::move(muon));

  std::cout << "rarexsec selectors chained successfully." << std::endl;
}
