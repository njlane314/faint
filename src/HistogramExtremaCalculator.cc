#include "faint/HistogramExtremaCalculator.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

#include "TH1.h"

namespace faint {
namespace plot {

void HistogramExtremaCalculator::set_histogram(const TH1& histogram) {
  histogram_ = &histogram;
}

double HistogramExtremaCalculator::maximum() const {
  if (!histogram_) {
    throw std::runtime_error(
        "HistogramExtremaCalculator::maximum histogram has not been set");
  }

  double max_value = -std::numeric_limits<double>::infinity();
  for (int bin = 1; bin <= histogram_->GetNbinsX(); ++bin) {
    max_value = std::max(max_value, histogram_->GetBinContent(bin));
  }

  return max_value;
}

double HistogramExtremaCalculator::maximum_with_error() const {
  if (!histogram_) {
    throw std::runtime_error(
        "HistogramExtremaCalculator::maximum_with_error histogram has not been set");
  }

  double max_value = -std::numeric_limits<double>::infinity();
  for (int bin = 1; bin <= histogram_->GetNbinsX(); ++bin) {
    max_value = std::max(max_value, histogram_->GetBinContent(bin) +
                                         histogram_->GetBinError(bin));
  }

  return max_value;
}

}  // namespace plot
}  // namespace faint
