#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TSystem.h"
#include "ROOT/RDataFrame.hxx"

#include <faint/Campaign.h>

void example_macro(const char* run_config = nullptr)
{
  if (gSystem->Load("libfaint_root")) {
    throw std::runtime_error("Failed to load libfaint_root library");
  }

  std::string config_path = run_config ? run_config : "";
  if (config_path.empty()) {
    config_path = faint::campaign::run_config_path();
  }

  ROOT::EnableImplicitMT();

  faint::campaign::Options options;
  options.beam = "numi-fhc";
  options.periods = {"run1"};
  options.ntuple_dir = faint::campaign::ntuple_directory(config_path);

  auto campaign = faint::campaign::Campaign::open(config_path, options);

  std::cout << "Loaded beam " << campaign.beam() << " for";
  for (const auto& p : campaign.periods()) {
    std::cout << ' ' << p;
  }
  std::cout << " with " << campaign.sample_keys().size() << " samples." << std::endl;

  for (const auto& key : campaign.sample_keys()) {
    auto final_count = campaign.final(key).Count();
    auto final_count_value = final_count.GetValue();
    std::cout << "Final selection entries for " << key << ": " << final_count_value << std::endl;
  }

  std::cout << "Total POT: " << campaign.pot() << std::endl;
  std::cout << "Total triggers: " << campaign.triggers() << std::endl;
}
