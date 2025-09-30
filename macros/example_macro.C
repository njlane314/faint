#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>
#include <faint/Dataset.h>
#include <faint/Log.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void example_macro() {
  try {
    ROOT::EnableImplicitMT();

    if (gSystem->Load("libfaint_root")) {
      throw std::runtime_error("Failed to load libfaint_root library");
    }

    const std::string config_path = "data/samples.json";

    faint::dataset::Options options;
    options.beam = "numi-fhc";
    options.periods = {"run1"};
    options.ntuple_dir = faint::dataset::ntuple_directory();

    auto dataset = faint::dataset::Dataset::open(config_path, options);

    std::cout << "Loaded beam " << dataset.beam() << " for";
    for (const auto& period : dataset.periods()) {
      std::cout << ' ' << period;
    }
    std::cout << " with " << dataset.sample_keys().size() << " samples." << std::endl;

    for (const auto& key : dataset.sample_keys()) {
      auto final_count = dataset.final(key).Count();
      std::cout << "Final selection entries for " << key << ": " << final_count.GetValue()
                << std::endl;
    }

    std::cout << "Total POT: " << dataset.pot() << std::endl;
    std::cout << "Total triggers: " << dataset.triggers() << std::endl;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
  }
}
