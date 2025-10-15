#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>

#include <rarexsec/proc/Env.h>
#include <rarexsec/Hub.h>
#include <rarexsec/proc/Snapshot.h>

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

        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();
        const auto samples = hub.simulation_entries(env.beamline, env.periods);

        rarexsec::snapshot::Options opt;
        opt.outdir = "snapshots";
        opt.tree = env.tree;
        std::string outfile = env.beamline;
        for (const auto& period : env.periods) {
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
