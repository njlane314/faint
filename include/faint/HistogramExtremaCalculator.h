#ifndef FAINT_HISTOGRAM_EXTREMA_CALCULATOR_H
#define FAINT_HISTOGRAM_EXTREMA_CALCULATOR_H

class TH1;

namespace faint {
namespace plot {

class HistogramExtremaCalculator {
 public:
  HistogramExtremaCalculator() = default;

  void set_histogram(const TH1& histogram);

  [[nodiscard]] double maximum() const;
  [[nodiscard]] double maximum_with_error() const;

 private:
  const TH1* histogram_{nullptr};
};

}  // namespace plot
}  // namespace faint

#endif  // FAINT_HISTOGRAM_EXTREMA_CALCULATOR_H
