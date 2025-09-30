#ifndef FAINT_MATRIX_HISTOGRAM_BUILDER_H
#define FAINT_MATRIX_HISTOGRAM_BUILDER_H

#include <memory>
#include <string>

class TH2D;
class TMatrixD;

namespace faint {
namespace plot {

class MatrixHistogramBuilder {
 public:
  MatrixHistogramBuilder() = default;

  void set_matrix(const TMatrixD& matrix);
  void set_template(const TH2D& histogram);
  void set_name(std::string name);

  [[nodiscard]] std::unique_ptr<TH2D> build() const;

 private:
  const TMatrixD* matrix_{nullptr};
  const TH2D* template_histogram_{nullptr};
  std::string name_{};
};

}  // namespace plot
}  // namespace faint

#endif  // FAINT_MATRIX_HISTOGRAM_BUILDER_H
