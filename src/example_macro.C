#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TSystem.h"

#include <faint/Campaign.h>

void example_macro()
{
  if (gSystem->Load("libfaint_root")) {
    throw std::runtime_error("Failed to load libfaint_root library");
  }

  faint::campaign::Options options;
  options.beam = "numi-fhc";
  options.periods = {"run1"};
  options.ntuple_dir = faint::campaign::ntuple_directory();

  auto campaign = faint::campaign::Campaign::open(faint::campaign::run_config_path(), options);

  std::cout << "Loaded beam " << campaign.beam() << " for";
  for (const auto& p : campaign.periods()) {
    std::cout << ' ' << p;
  }
  std::cout << " with " << campaign.sample_keys().size() << " samples." << std::endl;

  for (const auto& key : campaign.sample_keys()) {
    auto final_count = campaign.final(key).Count();
    std::cout << "Final selection entries for " << key << ": " << final_count.GetValue() << std::endl;
  }

  std::cout << "Total POT: " << campaign.pot() << std::endl;
  std::cout << "Total triggers: " << campaign.triggers() << std::endl;
}
