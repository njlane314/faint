#include "faint/MatrixHistogramBuilder.h"

#include <stdexcept>
#include <utility>

#include "TH2D.h"
#include "TMatrixD.h"

namespace faint {
namespace plot {

void MatrixHistogramBuilder::set_matrix(const TMatrixD& matrix) {
  matrix_ = &matrix;
}

void MatrixHistogramBuilder::set_template(const TH2D& histogram) {
  template_histogram_ = &histogram;
}

void MatrixHistogramBuilder::set_name(std::string name) { name_ = std::move(name); }

std::unique_ptr<TH2D> MatrixHistogramBuilder::build() const {
  if (!matrix_ || !template_histogram_) {
    throw std::runtime_error(
        "MatrixHistogramBuilder::build requires both matrix and template histogram");
  }

  const std::string histogram_name = name_.empty() ?
                                         std::string("h_") + template_histogram_->GetName() :
                                         name_;
  auto histogram = std::unique_ptr<TH2D>(
      static_cast<TH2D*>(template_histogram_->Clone(histogram_name.c_str())));

  const int nbins = histogram->GetNbinsX();
  for (int i = 1; i <= nbins; ++i) {
    for (int j = 1; j <= nbins; ++j) {
      histogram->SetBinContent(i, j, (*matrix_)[i - 1][j - 1]);
    }
  }

  return histogram;
}

}  // namespace plot
}  // namespace faint
