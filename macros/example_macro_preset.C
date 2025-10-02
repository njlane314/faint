#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>
#include <rarexsec/Hub.hh>
#include <rarexsec/Selection.hh>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void example_macro_preset() {
    try {
        ROOT::EnableImplicitMT();
        if (gSystem->Load("librarexsec") < 0) {
            throw std::runtime_error("Failed to load librarexsec library");
        }

        const std::string config_path = "data/samples.json";
        const std::string beamline = "numi-fhc";
        const std::vector<std::string> periods = {"run1"};

        rarexsec::Hub hub(config_path);

        const auto preset = rarexsec::selection::Preset::Baseline;
        std::cout << "Using preset: " << rarexsec::selection::to_string(preset) << "\n";

        auto summary = rarexsec::selection::evaluate(hub, beamline, periods, preset, "w_nominal");
        rarexsec::selection::print(summary);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}
