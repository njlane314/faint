#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMatrixD.h>
#include <faint/Dataset.h>
#include <faint/Log.h>
#include <faint/Types.h>
#include <faint/plot/ErrorBandBuilder.h>
#include <faint/syst/systematics_rdf.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
// Utility: find the first Monte Carlo sample key in the dataset.
std::string find_mc_sample_key(const faint::dataset::Dataset& dataset) {
  std::vector<std::string> keys = dataset.sample_keys(faint::SampleOrigin::kMonteCarlo);
  if (keys.empty()) {
    throw std::runtime_error("Dataset does not contain a Monte Carlo sample");
  }
  return keys.front();
}

// Utility: convert a TMatrixD into a TH2D histogram.
TH2D matrix_to_histogram(const TMatrixD& matrix, const std::string& name) {
  const int nbins = matrix.GetNrows();
  TH2D hist(name.c_str(), "Systematic covariance;Bin;Bin", nbins, 0.5, nbins + 0.5,
            nbins, 0.5, nbins + 0.5);
  for (int i = 0; i < nbins; ++i) {
    for (int j = 0; j < nbins; ++j) {
      hist.SetBinContent(i + 1, j + 1, matrix[i][j]);
    }
  }
  return hist;
}

}  // namespace

void systematics_error_band_example() {
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

    const auto sample_key = find_mc_sample_key(dataset);
    auto final_df = dataset.final(sample_key);

    const std::string observable = "contained_fraction";
    ROOT::RDF::TH1DModel hist_model{"h_cv", "Contained fraction;Contained fraction;Events", 20,
                                    0.0, 1.0};

    auto cv_result = faint::syst::cv_histogram(final_df, hist_model, observable);
    const TH1D* cv_hist = cv_result.GetPtr();

    std::vector<ROOT::RDF::RResultPtr<TH1D>> universe_results_storage;
    TMatrixD total_covariance(cv_hist->GetNbinsX(), cv_hist->GetNbinsX());
    total_covariance.Zero();

    const auto systematics = faint::syst::systematic_list_from_variables();
    for (const auto& syst : systematics) {
      auto universe_results =
          faint::syst::universe_histograms(final_df, syst, hist_model, observable);
      if (universe_results.empty()) {
        continue;
      }

      std::vector<const TH1D*> universe_histograms;
      universe_histograms.reserve(universe_results.size());
      for (auto& result : universe_results) {
        universe_histograms.push_back(result.GetPtr());
        universe_results_storage.emplace_back(result);
      }

      auto covariance =
          faint::syst::covariance_matrix_from_histograms(cv_hist, universe_histograms, syst.kind);
      total_covariance += covariance;
    }

    auto covariance_hist = matrix_to_histogram(total_covariance, "h_total_covariance");

    faint::plot::ErrorBandBuilder builder;
    builder.add_component(*cv_hist);
    builder.set_covariance(covariance_hist);
    auto error_band = builder.build("_syst");

    error_band->SetFillColorAlpha(kAzure + 1, 0.35);
    error_band->SetLineColor(kAzure + 2);
    error_band->SetMarkerSize(0.0);

    auto nominal_clone = std::unique_ptr<TH1D>(
        static_cast<TH1D*>(cv_hist->Clone("h_nominal_clone")));
    nominal_clone->SetLineColor(kBlack);
    nominal_clone->SetLineWidth(2);

    TCanvas canvas("c_systematics", "Histogram with systematic error band", 800, 600);
    error_band->Draw("E2");
    nominal_clone->Draw("HIST SAME");

    TLegend legend(0.55, 0.7, 0.88, 0.88);
    legend.AddEntry(nominal_clone.get(), "Nominal prediction", "l");
    legend.AddEntry(error_band.get(), "Total syst. error", "f");
    legend.Draw();

    canvas.Update();
    canvas.SaveAs("systematics_error_band.png");

    std::cout << "Plotted systematic error band for sample: " << sample_key << std::endl;
    std::cout << "Output saved to systematics_error_band.png" << std::endl;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
  }
}
