#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>
#include <rarexsec/Hub.hh>
#include <rarexsec/Plotter.hh>
#include <stdexcept>

void make_plots() {
    ROOT::EnableImplicitMT();
    if (gSystem->Load("librarexsec") < 0)
        throw std::runtime_error("Failed to load librexsec");

    const std::string config = "data/samples.json";
    const std::string beamline = "NuMI FHC";
    const std::vector<std::string> periods = {"run1"};

    rarexsec::Hub hub(config);
    auto mc = hub.simulation_entries("numi-fhc", periods);
    auto data = hub.data_entries("numi-fhc", periods);

    rarexsec::plot::Plotter plotter({.out_dir = "plots/numi-fhc/run1",
                                     .image_format = "png",
                                     .show_ratio = true,
                                     .use_log_y = false,
                                     .annotate_numbers = true,
                                     .overlay_signal = true,
                                     .signal_channels = {15, 16},
                                     .y_min = 0.,
                                     .y_max = -1.,
                                     .cuts = {},
                                     .beamline = beamline,
                                     .periods = periods});

    rarexsec::plot::Histogram1DSpec h_len{
        .name = "trk_len",
        .title = ";L_{trk} [cm];Events",
        .expr = "ROOT::VecOps::Max(track_length)",
        .weight = "w_nominal",
        .nbins = 50,
        .xmin = 0.,
        .xmax = 200.,
        .sel = rarexsec::selection::Preset::InclusiveMuCC};

    plotter.draw_stack_by_channel(h_len, mc, data);
}
