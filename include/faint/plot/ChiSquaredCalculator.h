#ifndef FAINT_CHI_SQUARED_CALCULATOR_H
#define FAINT_CHI_SQUARED_CALCULATOR_H

#include <optional>
#include <utility>
#include <vector>

#include "TMatrixDSym.h"

class TH1D;

namespace faint {
namespace plot {

class ChiSquaredCalculator {
 public:
  ChiSquaredCalculator() = default;

  void set_prediction(const TH1D& histogram);
  void set_data(const TH1D& histogram);

  void set_covariance(const TMatrixDSym& covariance);
  void clear_covariance();

  void set_skip_bins(std::vector<int> bins);
  void add_skip_bin(int bin);
  void clear_skip_bins();

  [[nodiscard]] std::pair<double, int> compute() const;

 private:
  const TH1D* prediction_{nullptr};
  const TH1D* data_{nullptr};
  std::optional<TMatrixDSym> covariance_{};
  std::vector<int> skip_bins_{};
};

}  // namespace plot
}  // namespace faint

#endif  // FAINT_CHI_SQUARED_CALCULATOR_H
