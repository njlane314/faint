#include "faint/plot/ErrorBandBuilder.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "TH1D.h"
#include "TH2D.h"

namespace faint {
namespace plot {

void ErrorBandBuilder::set_components(
    const std::vector<const TH1D*>& histograms) {
  components_ = histograms;
}

void ErrorBandBuilder::add_component(const TH1D& histogram) {
  components_.push_back(&histogram);
}

void ErrorBandBuilder::clear_components() { components_.clear(); }

void ErrorBandBuilder::set_covariance(const TH2D& covariance) {
  covariance_ = &covariance;
}

void ErrorBandBuilder::clear_covariance() { covariance_ = nullptr; }

std::unique_ptr<TH1D> ErrorBandBuilder::build(const std::string& name_suffix) const {
  if (components_.empty()) {
    throw std::runtime_error(
        "ErrorBandBuilder::build requires at least one component histogram");
  }

  const std::string clone_name =
      std::string(components_.front()->GetName()) + name_suffix;
  auto error_histogram = std::unique_ptr<TH1D>(
      static_cast<TH1D*>(components_.front()->Clone(clone_name.c_str())));
  error_histogram->Reset();

  for (int bin = 1; bin <= error_histogram->GetNbinsX(); ++bin) {
    double events = 0.0;
    double variance = 0.0;

    for (const auto* histogram : components_) {
      events += histogram->GetBinContent(bin);
      const double error = histogram->GetBinError(bin);
      variance += error * error;
    }

    if (covariance_) {
      variance += covariance_->GetBinContent(bin, bin);
    }

    error_histogram->SetBinContent(bin, events);
    error_histogram->SetBinError(bin, std::sqrt(variance));
  }

  return error_histogram;
}

}  // namespace plot
}  // namespace faint
