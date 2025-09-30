#include "faint/DataGraphBuilder.h"

#include <memory>

#include "TGraph.h"

namespace faint {
namespace plot {

void DataGraphBuilder::set_positions(const std::vector<double>& positions) {
  positions_ = positions;
}

void DataGraphBuilder::set_marker_style(int style) { marker_style_ = style; }

void DataGraphBuilder::set_marker_size(double size) { marker_size_ = size; }

std::unique_ptr<TGraph> DataGraphBuilder::build() const {
  auto graph = std::make_unique<TGraph>(positions_.size());

  for (int i = 0; i < graph->GetN(); ++i) {
    graph->SetPoint(i, positions_.at(i), 0.0);
  }

  graph->SetMarkerStyle(marker_style_);
  graph->SetMarkerSize(marker_size_);

  return graph;
}

}  // namespace plot
}  // namespace faint
