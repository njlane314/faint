#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>

#include <rarexsec/Hub.hh>
#include <rarexsec/Snapshot.hh>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void snapshot_macro() {
    try {
        ROOT::EnableImplicitMT();

        if (gSystem->Load("librarexsec") < 0) {
            throw std::runtime_error("Failed to load librexsec");
        }

        const std::string config_path = "data/samples.json";
        const std::string beamline = "numi-fhc";
        const std::vector<std::string> periods = {"run1"};

        rarexsec::Hub hub(config_path);

        rarexsec::snapshot::Options opt;
        opt.outdir = "snapshots";
        opt.tree = "analysis";

        auto outputs = rarexsec::snapshot::write(hub, beamline, periods, opt);

        if (outputs.empty()) {
            std::cout << "[snapshot] no files were written (no matching samples?).\n";
        } else {
            std::cout << "[snapshot] wrote " << outputs.size() << " file(s):\n";
            for (const auto& f : outputs) {
                std::cout << "  " << f << "\n";
            }
        }

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
    }
}
