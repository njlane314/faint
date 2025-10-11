#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>
#include <rarexsec/Hub.hh>
#include <rarexsec/Selection.hh>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void apply_inclusive_mucc_preset() {
    try {
        ROOT::EnableImplicitMT();
        if (gSystem->Load("librarexsec") < 0) {
            throw std::runtime_error("Failed to load librarexsec library");
        }

        const std::string config_path = "data/samples.json";
        const std::string beamline = "numi-fhc";
        const std::vector<std::string> periods = {"run1"};

        rarexsec::Hub hub(config_path);

        const auto preset = rarexsec::selection::Preset::InclusiveMuCC;
        auto preset_to_string = [](rarexsec::selection::Preset value) {
            switch (value) {
            case rarexsec::selection::Preset::Empty:
                return "Empty";
            case rarexsec::selection::Preset::Trigger:
                return "Trigger";
            case rarexsec::selection::Preset::Slice:
                return "Slice";
            case rarexsec::selection::Preset::Fiducial:
                return "Fiducial";
            case rarexsec::selection::Preset::Topology:
                return "Topology";
            case rarexsec::selection::Preset::Muon:
                return "Muon";
            case rarexsec::selection::Preset::InclusiveMuCC:
            default:
                return "InclusiveMuCC";
            }
        };

        std::cout << "Using preset: " << preset_to_string(preset) << "\n";

        const auto samples = hub.simulation_entries(beamline, periods);
        std::cout << "Found " << samples.size() << " simulation samples" << std::endl;

        for (const auto* entry : samples) {
            if (!entry) {
                continue;
            }

            auto node = rarexsec::selection::apply(entry->rnode(), preset, *entry);
            const auto selected = node.Count().GetValue();
            std::cout << "Sample '" << entry->file << "' selected entries: " << selected << std::endl;
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}
