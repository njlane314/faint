#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>
#include <faint/Campaign.h>
#include <faint/Log.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void example_macro() {
  try {
    faint::log::init();
    faint::log::set_level(faint::log::Level::kDebug);

    ROOT::EnableImplicitMT();

    if (gSystem->Load("libfaint_root")) {
      throw std::runtime_error("Failed to load libfaint_root library");
    }

    const std::string config_path = "data/samples.json";

    faint::campaign::Options options;
    options.beam = "numi-fhc";
    options.periods = {"run1"};
    options.ntuple_dir = faint::campaign::ntuple_directory();

    auto campaign = faint::campaign::Campaign::open(config_path, options);

    std::cout << "Loaded beam " << campaign.beam() << " for";
    for (const auto& period : campaign.periods()) {
      std::cout << ' ' << period;
    }
    std::cout << " with " << campaign.sample_keys().size() << " samples." << std::endl;

    for (const auto& key : campaign.sample_keys()) {
      auto final_count = campaign.final(key).Count();
      std::cout << "Final selection entries for " << key << ": " << final_count.GetValue()
                << std::endl;
    }

    std::cout << "Total POT: " << campaign.pot() << std::endl;
    std::cout << "Total triggers: " << campaign.triggers() << std::endl;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
  }
}
