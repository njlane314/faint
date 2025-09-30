#ifndef FAINT_DATA_GRAPH_BUILDER_H
#define FAINT_DATA_GRAPH_BUILDER_H

#include <memory>
#include <vector>

class TGraph;

namespace faint {
namespace plot {

class DataGraphBuilder {
 public:
  DataGraphBuilder() = default;

  void set_positions(const std::vector<double>& positions);
  void set_marker_style(int style);
  void set_marker_size(double size);

  [[nodiscard]] std::unique_ptr<TGraph> build() const;

 private:
  std::vector<double> positions_{};
  int marker_style_{23};
  double marker_size_{3.0};
};

}  // namespace plot
}  // namespace faint

#endif  // FAINT_DATA_GRAPH_BUILDER_H
