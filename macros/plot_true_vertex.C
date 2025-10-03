#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TSystem.h>
#include <TH1D.h>
#include <rarexsec/Hub.hh>

#include <array>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
struct VertexPlot {
    std::string column;
    std::string hist_name;
    std::string hist_title;
    std::string canvas_name;
    std::string output_path;
    int bins;
    double min;
    double max;
};
}

void plot_true_vertex() {
    ROOT::EnableImplicitMT();

    if (gSystem->Load("librarexsec") < 0) {
        throw std::runtime_error("Failed to load librexsec");
    }

    const std::string config_path = "data/samples.json";
    const std::string beamline = "numi-fhc";
    const std::vector<std::string> periods = {"run1"};

    rarexsec::Hub hub(config_path);
    const auto samples = hub.simulation_entries(beamline, periods);
    if (samples.empty()) {
        throw std::runtime_error("No simulation samples found for the requested configuration");
    }

    const std::string out_dir = "plots/" + beamline + "/" + periods.front();
    gSystem->mkdir(out_dir.c_str(), true);

    const std::array<VertexPlot, 3> plots = {
        VertexPlot{"neutrino_vertex_x", "h_true_vertex_x", "True neutrino vertex X;x^{true}_{#nu} [cm];Events", "c_true_vertex_x", out_dir + "/true_vertex_x.png", 50, 0., 256.},
        VertexPlot{"neutrino_vertex_y", "h_true_vertex_y", "True neutrino vertex Y;y^{true}_{#nu} [cm];Events", "c_true_vertex_y", out_dir + "/true_vertex_y.png", 50, -116., 116.},
        VertexPlot{"neutrino_vertex_z", "h_true_vertex_z", "True neutrino vertex Z;z^{true}_{#nu} [cm];Events", "c_true_vertex_z", out_dir + "/true_vertex_z.png", 80, 0., 1036.}
    };

    for (const auto& plot : plots) {
        TH1D total_hist(plot.hist_name.c_str(), plot.hist_title.c_str(), plot.bins, plot.min, plot.max);
        total_hist.Sumw2();

        std::size_t sample_index = 0;
        for (const auto* entry : samples) {
            if (!entry) {
                continue;
            }
            auto hist = entry->rnode().Histo1D({plot.hist_name + "_" + std::to_string(sample_index++), plot.hist_title.c_str(), plot.bins, plot.min, plot.max}, plot.column);
            auto* hist_ptr = hist->GetPtr();
            hist_ptr->SetDirectory(nullptr);
            total_hist.Add(hist_ptr);
        }

        if (total_hist.GetEntries() <= 0.) {
            throw std::runtime_error("Histogram " + plot.hist_name + " has no entries");
        }

        auto canvas = std::make_unique<TCanvas>(plot.canvas_name.c_str(), plot.canvas_name.c_str(), 800, 600);
        total_hist.SetLineWidth(2);
        total_hist.Draw("hist");
        canvas->SaveAs(plot.output_path.c_str());
    }
}
