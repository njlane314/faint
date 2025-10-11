#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <exception>
#include <iostream>

#include <rarexsec/Hub.hh>
#include <rarexsec/proc/Env.hh>

void plot_topology_variables() {
    try {
        ROOT::EnableThreadSafety();
        ROOT::EnableImplicitMT();

        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();
        const auto samples = hub.simulation_entries(env.beamline, env.periods);

        TH1D h_contained("h_contained_fraction",
                         "Contained Fraction;Contained Fraction;Weighted Events",
                         50, 0.0, 1.0);
        h_contained.Sumw2();
        h_contained.SetDirectory(nullptr);

        TH1D h_cluster("h_slice_cluster_fraction",
                       "Slice Cluster Fraction;Slice Cluster Fraction;Weighted Events",
                       50, 0.0, 1.0);
        h_cluster.Sumw2();
        h_cluster.SetDirectory(nullptr);

        for (const auto* entry : samples) {
            if (!entry) {
                continue;
            }

            auto node = entry->nominal.rnode();
            auto contained = node.Histo1D({"tmp_contained", "", 50, 0.0, 1.0},
                                          "contained_fraction",
                                          "w_nominal");
            auto cluster = node.Histo1D({"tmp_cluster", "", 50, 0.0, 1.0},
                                        "slice_cluster_fraction",
                                        "w_nominal");

            h_contained.Add(contained.GetPtr());
            h_cluster.Add(cluster.GetPtr());
        }

        TCanvas* canvas = new TCanvas("c_selection_variables",
                                      "Topology Selection Variables",
                                      800, 600);
        canvas->Divide(1, 2);

        canvas->cd(1);
        h_contained.Draw();

        canvas->cd(2);
        h_cluster.Draw();

        canvas->Update();

        ROOT::DisableImplicitMT();
    } catch (const std::exception& ex) {
        ROOT::DisableImplicitMT();
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}
