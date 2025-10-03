#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>
#include <ROOT/RDFHelpers.hxx>
#include <rarexsec/Hub.hh>
#include <rarexsec/Plotter.hh>

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

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

    rarexsec::plot::Options plot_options;
    plot_options.out_dir = out_dir;
    plot_options.image_format = "png";
    plot_options.show_ratio = false;
    plot_options.beamline = beamline;
    plot_options.periods = periods;

    rarexsec::plot::Plotter plotter(plot_options);

    const std::array<rarexsec::plot::H1Spec, 3> plots = {
        rarexsec::plot::H1Spec{
            .id = "true_vertex_x",
            .title = "True neutrino vertex X;x^{true}_{#nu} [cm];Events",
            .expr = "neutrino_vertex_x",
            .weight = "w_nominal",
            .nbins = 50,
            .xmin = 0.,
            .xmax = 256.,
            .sel = rarexsec::selection::Preset::Empty
        },
        rarexsec::plot::H1Spec{
            .id = "true_vertex_y",
            .title = "True neutrino vertex Y;y^{true}_{#nu} [cm];Events",
            .expr = "neutrino_vertex_y",
            .weight = "w_nominal",
            .nbins = 50,
            .xmin = -116.,
            .xmax = 116.,
            .sel = rarexsec::selection::Preset::Empty
        },
        rarexsec::plot::H1Spec{
            .id = "true_vertex_z",
            .title = "True neutrino vertex Z;z^{true}_{#nu} [cm];Events",
            .expr = "neutrino_vertex_z",
            .weight = "w_nominal",
            .nbins = 80,
            .xmin = 0.,
            .xmax = 1036.,
            .sel = rarexsec::selection::Preset::Empty
        }
    };

    for (const auto& plot_spec : plots) {
        plotter.draw_stack_by_channel(plot_spec, samples);
    }
}
