#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>
#include <rarexsec/Dataset.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void example_macro() {
  try {
    ROOT::EnableImplicitMT();

    if (gSystem->Load("librarexsec_root")) {
      throw std::runtime_error("Failed to load librarexsec_root library");
    }

    const std::string config_path = "data/samples.json";

    rarexsec::dataset::Options options;
    options.beam = "numi-fhc";
    options.periods = {"run1"};
    auto dataset = rarexsec::dataset::Dataset::open(config_path, options);

    std::cout << "Loaded beam " << dataset.beam() << " for";
    for (const auto& period : dataset.periods()) {
      std::cout << ' ' << period;
    }
    std::cout << " with " << dataset.sample_keys().size() << " samples." << std::endl;

    for (const auto& key : dataset.sample_keys()) {
      auto final_count = dataset.final(key).Count();
      auto eval_final_count = final_count.GetValue();
      std::cout << "Final selection entries for " << key << ": " << eval_final_count << std::endl;
    }

    std::cout << "Total POT: " << dataset.pot() << std::endl;
    std::cout << "Total triggers: " << dataset.triggers() << std::endl;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
  }
}
