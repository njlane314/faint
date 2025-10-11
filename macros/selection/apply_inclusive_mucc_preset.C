#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>

#include <rarexsec/Hub.hh>
#include <rarexsec/proc/DataModel.hh>
#include <rarexsec/proc/Env.hh>
#include <rarexsec/proc/Selection.hh>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void apply_inclusive_mucc_preset() {
    try {
        ROOT::EnableThreadSafety();
        ROOT::EnableImplicitMT();
    
        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();

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

        const auto samples = hub.simulation_entries(env.beamline, env.periods);
        std::cout << "Found " << samples.size() << " simulation samples" << std::endl;

        for (const auto* entry : samples) {
            if (!entry) {
                continue;
            }

            auto node = rarexsec::selection::apply(entry->rnode(), preset, *entry);
            const auto selected = node.Count().GetValue();
            std::cout << "Sample '" << entry->file << "' selected entries: " << selected << std::endl;
        }

        const auto eval = rarexsec::selection::evaluate(
            samples,
            [](int ch) {
                const auto channel = static_cast<rarexsec::Channel>(ch);
                switch (channel) {
                case rarexsec::Channel::MuCC0pi_ge1p:
                case rarexsec::Channel::MuCC1pi:
                case rarexsec::Channel::MuCCPi0OrGamma:
                case rarexsec::Channel::MuCCNpi:
                case rarexsec::Channel::MuCCOther:
                    return true;
                default:
                    return false;
                }
            },
            preset);

        std::cout << "Selection evaluation:\n"
                  << "  Denominator (signal truth): " << eval.denom << '\n'
                  << "  Selected (all): " << eval.selected << '\n'
                  << "  Selected signal: " << eval.numer << '\n'
                  << "  Efficiency: " << eval.efficiency() << '\n'
                  << "  Purity: " << eval.purity() << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}
