#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>
#include <TH1D.h>

#include <faint/Dataset.h>
#include <faint/Log.h>
#include <faint/Samples.h>
#include <faint/plot/StackedHistogram.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
ROOT::RDF::TH1DModel make_hist_model(const std::string& name,
                                     const std::string& title) {
  return ROOT::RDF::TH1DModel(name.c_str(), title.c_str(), 40, 0.0, 400.0);
}

std::vector<Color_t> default_palette() {
  return {kAzure + 1, kOrange + 7, kSpring + 5, kViolet + 5, kGray + 2};
}

std::string prettify_label(const std::string& key) {
  std::string label = key;
  std::replace(label.begin(), label.end(), '_', ' ');
  return label;
}
}  // namespace

void stacked_histogram_example() {
  try {
    ROOT::EnableImplicitMT();

    if (gSystem->Load("libfaint_root")) {
      throw std::runtime_error("Failed to load libfaint_root library");
    }

    const std::string config_path = "data/samples.json";

    faint::dataset::Options options;
    options.beam = "numi-fhc";
    options.periods = {"run1"};
    options.ntuple_dir = faint::dataset::ntuple_directory();

    auto dataset = faint::dataset::Dataset::open(config_path, options);

    faint::plot::StackedHistogram plot("muon_track_length");
    plot.set_x_axis_title("Muon candidate track length [cm]");
    plot.set_y_axis_title("Events");
    plot.set_legend_header(dataset.beam());
    plot.set_legend_position(0.65, 0.55, 0.88, 0.88);
    plot.set_annotate_yields(true);

    const std::vector<Color_t> palette = default_palette();
    std::size_t color_index = 0;

    auto add_background_samples = [&](faint::SampleOrigin origin) {
      auto keys = dataset.sample_keys(origin);
      for (const auto& key : keys) {
        const std::string hist_name = "hist_" + key;
        auto node = dataset.final(key);
        auto hist = node.Histo1D(make_hist_model(hist_name, key), "track_length",
                                 faint::dataset::col::Weight);
        const Color_t color = palette[color_index % palette.size()];
        ++color_index;
        plot.add_background(*hist, prettify_label(key), color, 1001);
      }
    };

    add_background_samples(faint::SampleOrigin::kMonteCarlo);
    add_background_samples(faint::SampleOrigin::kExternal);
    add_background_samples(faint::SampleOrigin::kDirt);

    auto data_keys = dataset.sample_keys(faint::SampleOrigin::kData);
    if (!data_keys.empty()) {
      const auto& key = data_keys.front();
      const std::string hist_name = "data_" + key;
      auto node = dataset.final(key);
      auto hist =
          node.Histo1D(make_hist_model(hist_name, key), "track_length",
                       faint::dataset::col::Weight);
      plot.set_data(*hist, "Data", kBlack, 20);
    }

    plot.add_cut(200.0, faint::plot::StackedHistogram::CutDirection::kGreaterThan,
                 "Track length > 200 cm", kRed + 1);

    plot.draw_and_save("png");

    std::cout << "Stacked histogram saved to plots/muon_track_length.png" << std::endl;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
  }
}
