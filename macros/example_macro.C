#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>
#include <rarexsec/Hub.hh>

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
        const std::string beamline = "numi-fhc";
        const std::vector<std::string> periods = {"run1"};

        rarexsec::Hub hub(config_path);
        const auto samples = hub.simulation(beamline, periods);

        std::cout << "Loaded beamline " << beamline << " for";
        for (const auto& period : periods) {
            std::cout << ' ' << period;
        }
        std::cout << " with " << samples.size() << " simulation samples." << std::endl;

        auto origin_to_string = [](rarexsec::sample::origin kind) {
            switch (kind) {
                case rarexsec::sample::origin::data:
                    return "data";
                case rarexsec::sample::origin::beam:
                    return "beam";
                case rarexsec::sample::origin::strangeness:
                    return "strangeness";
                case rarexsec::sample::origin::ext:
                    return "ext";
                case rarexsec::sample::origin::dirt:
                    return "dirt";
                case rarexsec::sample::origin::unknown:
                default:
                    return "unknown";
            }
        };

        double total_pot = 0.0;
        double total_pot_eff = 0.0;
        double total_triggers = 0.0;
        double total_triggers_eff = 0.0;

        for (const auto* entry : samples) {
            if (!entry) {
                continue;
            }

            std::cout << "Sample kind '" << origin_to_string(entry->kind) << "' from file "
                      << entry->file << std::endl;

            auto final_count = entry->rnode().Count();
            auto eval_final_count = final_count.GetValue();
            std::cout << "  Final selection entries: " << eval_final_count << std::endl;

            for (const auto& detvar : entry->detvars) {
                if (!detvar.second.node) {
                    continue;
                }
                auto detvar_count = detvar.second.rnode().Count().GetValue();
                std::cout << "  Detector variation '" << detvar.first
                          << "' entries: " << detvar_count << std::endl;
            }

            total_pot += entry->pot;
            total_pot_eff += (entry->pot_eff > 0.0) ? entry->pot_eff : entry->pot;
            total_triggers += entry->trig;
            total_triggers_eff += (entry->trig_eff > 0.0) ? entry->trig_eff : entry->trig;
        }

        std::cout << "Total POT (nominal): " << total_pot << std::endl;
        std::cout << "Total POT (effective): " << total_pot_eff << std::endl;
        std::cout << "Total triggers (nominal): " << total_triggers << std::endl;
        std::cout << "Total triggers (effective): " << total_triggers_eff << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}
