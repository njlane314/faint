#include "faint/plot/ChiSquaredCalculator.h"

#include <algorithm>
#include <stdexcept>

#include "TMatrixDSym.h"
#include "TH1D.h"

namespace faint {
namespace plot {
namespace {

bool bin_should_be_skipped(int bin, const std::vector<int>& skip_bins) {
  return std::find(skip_bins.begin(), skip_bins.end(), bin) != skip_bins.end();
}

}  // namespace

void ChiSquaredCalculator::set_prediction(const TH1D& histogram) {
  prediction_ = &histogram;
}

void ChiSquaredCalculator::set_data(const TH1D& histogram) { data_ = &histogram; }

void ChiSquaredCalculator::set_covariance(const TMatrixDSym& covariance) {
  covariance_ = covariance;
}

void ChiSquaredCalculator::clear_covariance() { covariance_.reset(); }

void ChiSquaredCalculator::set_skip_bins(std::vector<int> bins) {
  skip_bins_ = std::move(bins);
}

void ChiSquaredCalculator::add_skip_bin(int bin) { skip_bins_.push_back(bin); }

void ChiSquaredCalculator::clear_skip_bins() { skip_bins_.clear(); }

std::pair<double, int> ChiSquaredCalculator::compute() const {
  if (!prediction_ || !data_) {
    throw std::runtime_error(
        "ChiSquaredCalculator::compute requires both prediction and data histograms");
  }

  const int nbins = data_->GetNbinsX();
  std::vector<int> nonzero_bins;
  nonzero_bins.reserve(nbins);

  for (int bin = 1; bin <= nbins; ++bin) {
    if (prediction_->GetBinContent(bin) > 0.0 &&
        !bin_should_be_skipped(bin, skip_bins_)) {
      nonzero_bins.push_back(bin);
    }
  }

  if (nonzero_bins.empty()) {
    return {0.0, 0};
  }

  TMatrixDSym covariance_matrix(nonzero_bins.size());

  if (covariance_) {
    if (covariance_->GetNcols() != nbins) {
      throw std::invalid_argument(
          "ChiSquaredCalculator::compute covariance matrix dimensions do not match "
          "histogram bins");
    }

    for (std::size_t i = 0; i < nonzero_bins.size(); ++i) {
      for (std::size_t j = 0; j < nonzero_bins.size(); ++j) {
        covariance_matrix[i][j] =
            (*covariance_)[nonzero_bins.at(i) - 1][nonzero_bins.at(j) - 1];
      }
    }
  }

  TMatrixDSym statistical_covariance(nonzero_bins.size());
  for (std::size_t i = 0; i < nonzero_bins.size(); ++i) {
    const int bin = nonzero_bins.at(i);
    const double prediction_error = prediction_->GetBinError(bin);
    const double data_error = data_->GetBinError(bin);
    statistical_covariance[i][i] = prediction_error * prediction_error +
                                   data_error * data_error;
  }

  covariance_matrix += statistical_covariance;
  covariance_matrix.Invert();

  double chi2 = 0.0;
  for (std::size_t i = 0; i < nonzero_bins.size(); ++i) {
    const int bin_i = nonzero_bins.at(i);
    const double diff_i =
        prediction_->GetBinContent(bin_i) - data_->GetBinContent(bin_i);
    for (std::size_t j = 0; j < nonzero_bins.size(); ++j) {
      const int bin_j = nonzero_bins.at(j);
      const double diff_j =
          prediction_->GetBinContent(bin_j) - data_->GetBinContent(bin_j);
      chi2 += diff_i * covariance_matrix[i][j] * diff_j;
    }
  }

  return {chi2, static_cast<int>(nonzero_bins.size())};
}

}  // namespace plot
}  // namespace faint
