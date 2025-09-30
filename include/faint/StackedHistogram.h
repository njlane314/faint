#ifndef FAINT_STACKED_HISTOGRAM_H
#define FAINT_STACKED_HISTOGRAM_H

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "Rtypes.h"
#include "TAttLine.h"
#include "TColor.h"

class TCanvas;
class TH1;
class TObject;

namespace faint {
namespace plot {

class StackedHistogram {
 public:
  enum class CutDirection { kLessThan, kGreaterThan };

  struct Cut {
    double threshold{0.0};
    CutDirection direction{CutDirection::kGreaterThan};
    std::string label{};
    Color_t color{kRed};
  };

  StackedHistogram(std::string plot_name, std::string output_directory = "plots");
  ~StackedHistogram();

  void set_x_axis_title(std::string title);
  void set_y_axis_title(std::string title);
  void set_log_y(bool value);
  void set_y_axis_range(double minimum, double maximum);
  void reset_y_axis_range();

  void set_legend_position(double x1, double y1, double x2, double y2);
  void set_legend_columns(int columns);
  void set_legend_text_size(double size);
  void set_legend_header(std::optional<std::string> header);

  void set_annotate_yields(bool value);

  void add_background(const TH1& hist, std::string label, Color_t color,
                     Style_t fill_style = 1001);
  void clear_backgrounds();

  void set_data(const TH1& hist, std::string label = "Data",
               Color_t color = kBlack, Style_t marker_style = 20);
  void clear_data();

  void set_signal(const TH1& hist, std::string label,
                 Color_t color = kGreen + 2, Style_t line_style = kDashed,
                 double scale = 1.0, int line_width = 2);
  void clear_signal();

  void add_cut(double threshold, CutDirection direction,
              std::string label = std::string(), Color_t color = kRed);
  void clear_cuts();

  void draw(TCanvas& canvas);
  void draw_and_save(const std::string& format = "png");

 private:
  struct BackgroundComponent {
    std::string label;
    std::unique_ptr<TH1> histogram;
    Color_t color{kGray + 1};
    Style_t fill_style{1001};
  };

  struct DataComponent {
    std::string label;
    std::unique_ptr<TH1> histogram;
    Color_t color{kBlack};
    Style_t marker_style{20};
  };

  struct SignalComponent {
    std::string label;
    std::unique_ptr<TH1> histogram;
    Color_t color{kGreen + 2};
    Style_t line_style{kDashed};
    double scale{1.0};
    int line_width{2};
  };

  std::unique_ptr<TH1> clone_histogram(const TH1& hist,
                                       const std::string& suffix) const;
  std::string format_yield(double value, int precision = 1) const;
  void set_global_style() const;
  void draw_cuts(double max_y, double x_min, double x_max);

  std::string plot_name_;
  std::string output_directory_;
  std::string x_axis_title_;
  std::string y_axis_title_{"Events"};
  bool use_log_y_{false};
  bool has_y_range_{false};
  double y_min_{0.0};
  double y_max_{0.0};
  double legend_x1_{0.62};
  double legend_y1_{0.6};
  double legend_x2_{0.88};
  double legend_y2_{0.88};
  int legend_columns_{1};
  double legend_text_size_{0.04};
  std::optional<std::string> legend_header_{};
  bool annotate_yields_{true};

  std::vector<BackgroundComponent> backgrounds_;
  std::optional<DataComponent> data_;
  std::optional<SignalComponent> signal_;
  std::vector<Cut> cuts_;
  std::vector<std::unique_ptr<TObject>> overlays_;
};

}  // namespace plot
}  // namespace faint

#endif  // FAINT_STACKED_HISTOGRAM_H
