#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>

#include <rarexsec/Hub.hh>
#include <rarexsec/Snapshot.hh>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void write_simulation_snapshots() {
    try {
        ROOT::EnableImplicitMT();

        if (gSystem->Load("librarexsec") < 0) {
            throw std::runtime_error("Failed to load librexsec");
        }

        const std::string config_path = "data/samples.json";
        const std::string beamline = "numi-fhc";
        const std::vector<std::string> periods = {"run1"};

        rarexsec::Hub hub(config_path);
        const auto samples = hub.simulation_entries(beamline, periods);

        rarexsec::snapshot::Options opt;
        opt.outdir = "snapshots";
        opt.tree = "analysis";
        std::string outfile = beamline;
        for (const auto& period : periods) {
            outfile += "_" + period;
        }
        outfile += ".root";
        opt.outfile = outfile;

        auto outputs = rarexsec::snapshot::write(samples, opt);

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
