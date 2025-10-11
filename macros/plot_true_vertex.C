#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>
#include <rarexsec/Hub.hh>
#include <rarexsec/Plotter.hh>

#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void plot_true_vertex() {
    std::cout << "Starting plot_true_vertex macro" << std::endl;
    ROOT::EnableImplicitMT();

    std::cout << "Loading librexsec library" << std::endl;
    if (gSystem->Load("librarexsec") < 0) {
        throw std::runtime_error("Failed to load librexsec");
    }

    const std::string config_path = "data/samples.json";
    const std::string beamline = "numi-fhc";
    const std::vector<std::string> periods = {"run1"};

    rarexsec::Hub hub(config_path);
    std::cout << "Retrieving simulation entries for beamline '" << beamline << "'" << std::endl;
    const auto samples = hub.simulation_entries(beamline, periods);
    if (samples.empty()) {
        throw std::runtime_error("No simulation samples found for the requested configuration");
    }
    std::cout << "Retrieved " << samples.size() << " simulation samples" << std::endl;

    const std::string out_dir = "plots/" + beamline + "/" + periods.front();
    std::cout << "Creating output directory: " << out_dir << std::endl;
    gSystem->mkdir(out_dir.c_str(), true);

    rarexsec::plot::Options plot_options;
    plot_options.out_dir = out_dir;
    plot_options.image_format = "pdf";
    plot_options.show_ratio = false;
    plot_options.overlay_signal = false;
    plot_options.legend_on_top = true;
    plot_options.legend_split = 0.85;
    plot_options.y_title = "Events";
    plot_options.beamline = beamline;
    plot_options.periods = periods;

    rarexsec::plot::Plotter plotter(plot_options);

    const std::array<rarexsec::plot::Histogram1DSpec, 3> plots = {
        rarexsec::plot::Histogram1DSpec{
            .id = "true_vertex_x",
            .title = ";x^{true}_{#nu} [cm];Events",
            .expr = "neutrino_vertex_x",
            .weight = "w_nominal",
            .nbins = 50,
            .xmin = 0.,
            .xmax = 256.,
            .sel = rarexsec::selection::Preset::Empty},
        rarexsec::plot::Histogram1DSpec{
            .id = "true_vertex_y",
            .title = ";y^{true}_{#nu} [cm];Events",
            .expr = "neutrino_vertex_y",
            .weight = "w_nominal",
            .nbins = 50,
            .xmin = -116.,
            .xmax = 116.,
            .sel = rarexsec::selection::Preset::Empty},
        rarexsec::plot::Histogram1DSpec{
            .id = "true_vertex_z",
            .title = ";z^{true}_{#nu} [cm];Events",
            .expr = "neutrino_vertex_z",
            .weight = "w_nominal",
            .nbins = 80,
            .xmin = 0.,
            .xmax = 1036.,
            .sel = rarexsec::selection::Preset::Empty}};

    for (const auto& plot_spec : plots) {
        std::cout << "Drawing plot: " << plot_spec.id << std::endl;
        plotter.draw_stack_by_channel(plot_spec, samples);
    }

    std::cout << "Finished plot_true_vertex macro" << std::endl;
}
