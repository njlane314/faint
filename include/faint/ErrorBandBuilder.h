#ifndef FAINT_ERROR_BAND_BUILDER_H
#define FAINT_ERROR_BAND_BUILDER_H

#include <memory>
#include <vector>

class TH1D;
class TH2D;

namespace faint {
namespace plot {

class ErrorBandBuilder {
 public:
  ErrorBandBuilder() = default;

  void set_components(const std::vector<const TH1D*>& histograms);
  void add_component(const TH1D& histogram);
  void clear_components();

  void set_covariance(const TH2D& covariance);
  void clear_covariance();

  [[nodiscard]] std::unique_ptr<TH1D> build(const std::string& name_suffix =
                                                "_errors") const;

 private:
  std::vector<const TH1D*> components_{};
  const TH2D* covariance_{nullptr};
};

}  // namespace plot
}  // namespace faint

#endif  // FAINT_ERROR_BAND_BUILDER_H
