#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TSystem.h"

#include <faint/study.h>

void example_macro()
{
  if (gSystem->Load("libfaint_root")) {
    throw std::runtime_error("Failed to load libfaint_root library");
  }

  faint::study::Options options;
  options.beam = "numi-fhc";
  options.periods = {"run1"};
  options.ntuple_dir = faint::study::ntuple_directory();

  auto study = faint::study::Study::open(faint::study::run_config_path(), options);

  std::cout << "Loaded beam " << study.beam() << " for";
  for (const auto& p : study.periods()) {
    std::cout << ' ' << p;
  }
  std::cout << " with " << study.sample_keys().size() << " samples." << std::endl;

  for (const auto& key : study.sample_keys()) {
    auto final_count = study.final(key).Count();
    std::cout << "Final selection entries for " << key << ": " << final_count.GetValue() << std::endl;
  }

  std::cout << "Total POT: " << study.pot() << std::endl;
  std::cout << "Total triggers: " << study.triggers() << std::endl;
}
