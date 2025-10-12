#include <ROOT/RDataFrame.hxx>
#include <exception>
#include <iostream>
#include <vector>

#include <rarexsec/Hub.hh>
#include <rarexsec/proc/Env.hh>
#include <rarexsec/plot/Descriptors.hh>
#include <rarexsec/plot/Plotter.hh>

void plot_topology_variables() {
    try {
        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();
        const auto mc_samples = hub.simulation_entries(env.beamline, env.periods);

        std::cout << "Loaded beamline " << env.beamline << " for";
        for (const auto& period : env.periods) {
            std::cout << ' ' << period;
        }
        std::cout << " with " << mc_samples.size() << " simulation samples." << std::endl;

        rarexsec::plot::Options opt;
        opt.out_dir = "plots/selection";
        opt.use_log_y = true;
        opt.overlay_signal = true;
        opt.image_format = "png";
        opt.legend_on_top = true;
        opt.beamline = env.beamline;
        opt.periods = env.periods;
        opt.analysis_region_label = "Empty Selection";

        rarexsec::plot::Plotter plotter(opt);

        rarexsec::plot::Histogram1DSpec contained;
        contained.id = "contained_fraction";
        contained.title = ";Contained Fraction;Events";
        contained.nbins = 100;
        contained.xmin = 0.0;
        contained.xmax = 1.0;
        contained.sel = rarexsec::selection::Preset::Empty;

        rarexsec::plot::Histogram1DSpec cluster = contained;
        cluster.id = "slice_cluster_fraction";
        cluster.title = ";Slice Cluster Fraction;Events";

        const std::vector<rarexsec::plot::Histogram1DSpec> specs = {contained, cluster};
        for (const auto& spec : specs) {
            plotter.draw_stack_by_channel(spec, mc_samples);
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}
