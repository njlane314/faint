#include "faint/StackedHistogram.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "TArrow.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "THStack.h"

#include "faint/Log.h"
#include "faint/PlotStyle.h"

namespace faint {
namespace plot {

namespace {
constexpr double kDefaultCanvasWidth = 800;
constexpr double kDefaultCanvasHeight = 600;
}  // namespace

StackedHistogram::StackedHistogram(std::string plot_name,
                                   std::string output_directory)
    : plot_name_(std::move(plot_name)),
      output_directory_(std::move(output_directory)) {}

StackedHistogram::~StackedHistogram() = default;

void StackedHistogram::set_x_axis_title(std::string title) {
  x_axis_title_ = std::move(title);
}

void StackedHistogram::set_y_axis_title(std::string title) {
  y_axis_title_ = std::move(title);
}

void StackedHistogram::set_log_y(bool value) { use_log_y_ = value; }

void StackedHistogram::set_y_axis_range(double minimum, double maximum) {
  has_y_range_ = true;
  y_min_ = minimum;
  y_max_ = maximum;
}

void StackedHistogram::reset_y_axis_range() { has_y_range_ = false; }

void StackedHistogram::set_legend_position(double x1, double y1, double x2,
                                         double y2) {
  legend_x1_ = x1;
  legend_y1_ = y1;
  legend_x2_ = x2;
  legend_y2_ = y2;
}

void StackedHistogram::set_legend_columns(int columns) {
  legend_columns_ = std::max(1, columns);
}

void StackedHistogram::set_legend_text_size(double size) { legend_text_size_ = size; }

void StackedHistogram::set_legend_header(std::optional<std::string> header) {
  legend_header_ = std::move(header);
}

void StackedHistogram::set_annotate_yields(bool value) { annotate_yields_ = value; }

void StackedHistogram::add_background(const TH1& hist, std::string label,
                                     Color_t color, Style_t fill_style) {
  std::string suffix = label.empty() ? "background" : label;
  std::transform(suffix.begin(), suffix.end(), suffix.begin(), [](unsigned char ch) {
    return (std::isalnum(ch) || ch == '_') ? static_cast<char>(ch) : '_';
  });
  suffix += "_" + std::to_string(backgrounds_.size());

  backgrounds_.push_back(BackgroundComponent{std::move(label),
                                             clone_histogram(hist, suffix), color,
                                             fill_style});
}

void StackedHistogram::clear_backgrounds() { backgrounds_.clear(); }

void StackedHistogram::set_data(const TH1& hist, std::string label, Color_t color,
                               Style_t marker_style) {
  data_ = DataComponent{std::move(label), clone_histogram(hist, "data"), color,
                        marker_style};
}

void StackedHistogram::clear_data() { data_.reset(); }

void StackedHistogram::set_signal(const TH1& hist, std::string label,
                                 Color_t color, Style_t line_style,
                                 double scale, int line_width) {
  signal_ = SignalComponent{std::move(label), clone_histogram(hist, "signal"),
                            color, line_style, scale, line_width};
}

void StackedHistogram::clear_signal() { signal_.reset(); }

void StackedHistogram::add_cut(double threshold, CutDirection direction,
                              std::string label, Color_t color) {
  cuts_.push_back(Cut{threshold, direction, std::move(label), color});
}

void StackedHistogram::clear_cuts() { cuts_.clear(); }

void StackedHistogram::draw(TCanvas& canvas) {
  set_global_style();

  canvas.cd();
  canvas.Clear();
  canvas.SetLogy(use_log_y_);

  overlays_.clear();

  THStack stack((plot_name_ + "_stack").c_str(), plot_name_.c_str());

  std::vector<BackgroundComponent*> background_ordering;
  background_ordering.reserve(backgrounds_.size());
  for (auto& background : backgrounds_) {
    background_ordering.push_back(&background);
  }

  std::stable_sort(background_ordering.begin(), background_ordering.end(),
                   [](const BackgroundComponent* lhs,
                      const BackgroundComponent* rhs) {
                     return lhs->histogram->Integral() <
                            rhs->histogram->Integral();
                   });
  std::reverse(background_ordering.begin(), background_ordering.end());

  double max_y = 0.0;
  for (auto* component : background_ordering) {
    auto* hist = component->histogram.get();
    hist->SetDirectory(nullptr);
    hist->SetFillColor(component->color);
    hist->SetFillStyle(component->fill_style);
    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
    stack.Add(hist, "HIST");
    max_y = std::max(max_y, hist->GetMaximum());
  }

  bool has_stack = !background_ordering.empty();
  if (has_stack) {
    stack.Draw("HIST");
  } else if (data_) {
    data_->histogram->Draw("E1");
    max_y = std::max(max_y, data_->histogram->GetMaximum());
  } else if (signal_) {
    signal_->histogram->Draw("HIST");
    max_y = std::max(max_y, signal_->histogram->GetMaximum());
  } else {
    faint::log::warn("StackedHistogram::draw", "No histograms available to draw");
    canvas.Update();
    return;
  }

  TH1* frame = has_stack ? stack.GetHistogram() : nullptr;
  if (!frame) {
    frame = has_stack ? background_ordering.front()->histogram.get()
                      : data_ ? data_->histogram.get()
                              : signal_ ? signal_->histogram.get() : nullptr;
  }

  double background_yield_total = 0.0;
  for (auto* component : background_ordering) {
    background_yield_total += component->histogram->Integral(0, component->histogram->GetNbinsX() + 1);
  }

  if (frame) {
    frame->GetXaxis()->SetTitle(x_axis_title_.c_str());
    frame->GetYaxis()->SetTitle(y_axis_title_.c_str());
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitleOffset(1.2);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
  }

  double dynamic_max = max_y;

  if (signal_) {
    auto* hist = signal_->histogram.get();
    const double scale = signal_->scale;
    if (scale != 1.0) {
      hist->Scale(scale);
    }
    hist->SetLineColor(signal_->color);
    hist->SetLineStyle(signal_->line_style);
    hist->SetLineWidth(signal_->line_width);
    hist->SetFillStyle(0);
    hist->Draw("HIST SAME");
    dynamic_max = std::max(dynamic_max, hist->GetMaximum());
    if (scale != 1.0) {
      hist->Scale(1.0 / scale);
    }
  }

  if (data_) {
    auto* hist = data_->histogram.get();
    hist->SetLineColor(data_->color);
    hist->SetMarkerColor(data_->color);
    hist->SetMarkerStyle(data_->marker_style);
    hist->SetLineWidth(1);
    hist->Draw("E1 SAME");
    dynamic_max = std::max(dynamic_max, hist->GetMaximum());
  }

  double y_min = use_log_y_ ? 0.1 : 0.0;
  double y_max = dynamic_max > 0.0
                      ? dynamic_max * (use_log_y_ ? 10.0 : 1.25)
                      : (use_log_y_ ? 10.0 : 1.0);

  if (has_y_range_) {
    y_min = y_min_;
    y_max = y_max_;
  }

  if (frame) {
    frame->SetMinimum(y_min);
    frame->SetMaximum(y_max);
  } else if (data_) {
    data_->histogram->GetYaxis()->SetRangeUser(y_min, y_max);
  }

  auto legend = std::make_unique<TLegend>(legend_x1_, legend_y1_, legend_x2_, legend_y2_);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetNColumns(legend_columns_);
  legend->SetTextSize(legend_text_size_);
  if (legend_header_) {
    legend->SetHeader(legend_header_->c_str(), "C");
  }

  auto label_builder = [&](const std::string& base_label, double yield) {
    if (!annotate_yields_) {
      return base_label;
    }
    std::ostringstream ss;
    ss << base_label;
    if (!base_label.empty()) {
      ss << " ";
    }
    ss << format_yield(yield, 2);
    return ss.str();
  };

  for (auto* component : background_ordering) {
    legend->AddEntry(component->histogram.get(),
                     label_builder(component->label,
                                   component->histogram->Integral(0, component->histogram->GetNbinsX() + 1))
                         .c_str(),
                     "f");
  }

  if (signal_) {
    legend->AddEntry(signal_->histogram.get(),
                     label_builder(signal_->label,
                                   signal_->histogram->Integral(0, signal_->histogram->GetNbinsX() + 1))
                         .c_str(),
                     "l");
  }

  if (data_) {
    legend->AddEntry(data_->histogram.get(),
                     label_builder(data_->label,
                                   data_->histogram->Integral(0, data_->histogram->GetNbinsX() + 1))
                         .c_str(),
                     "lep");
  }

  legend->Draw();
  overlays_.push_back(std::unique_ptr<TObject>(legend.release()));

  if (frame) {
    draw_cuts(y_max, frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
  }

  if (has_stack) {
    stack.SetMaximum(y_max);
    stack.SetMinimum(y_min);
  }

  if (frame) {
    std::ostringstream title;
    title << plot_name_;
    if (annotate_yields_ && background_yield_total > 0.0) {
      title << " (Bkg: " << format_yield(background_yield_total, 2) << ")";
    }
    const std::string text = title.str();
    auto watermark = std::make_unique<TLatex>(canvas.GetLeftMargin() + 0.01,
                                              1.0 - canvas.GetTopMargin() + 0.01,
                                              text.c_str());
    watermark->SetNDC();
    watermark->SetTextFont(62);
    watermark->SetTextSize(0.045);
    watermark->SetTextAlign(13);
    watermark->Draw();
    overlays_.push_back(std::unique_ptr<TObject>(watermark.release()));
  }

  canvas.RedrawAxis();
  canvas.Update();
}

void StackedHistogram::draw_and_save(const std::string& format) {
  set_global_style();
  gROOT->SetBatch(kTRUE);
  if (gSystem) {
    gSystem->mkdir(output_directory_.c_str(), true);
  }

  TCanvas canvas(plot_name_.c_str(), plot_name_.c_str(),
                 static_cast<int>(kDefaultCanvasWidth),
                 static_cast<int>(kDefaultCanvasHeight));

  draw(canvas);

  const std::string out_path = output_directory_ + "/" + plot_name_ + "." + format;
  if (format == "pdf") {
    canvas.SaveAs(out_path.c_str());
  } else {
    std::unique_ptr<TImage> image(TImage::Create());
    if (!image) {
      faint::log::warn("StackedHistogram::draw_and_save", "Failed to create ROOT image object");
      return;
    }
    image->FromPad(&canvas);
    image->WriteImage(out_path.c_str());
  }
}

std::unique_ptr<TH1> StackedHistogram::clone_histogram(const TH1& hist,
                                                       const std::string& suffix) const {
  std::string base_name = hist.GetName();
  if (base_name.empty()) {
    base_name = plot_name_;
  }
  std::string clone_name = base_name + "_" + suffix;
  TObject* cloned = hist.Clone(clone_name.c_str());
  auto* cloned_hist = dynamic_cast<TH1*>(cloned);
  if (!cloned_hist) {
    throw std::runtime_error("Failed to clone histogram for stacked plot");
  }
  cloned_hist->SetDirectory(nullptr);
  return std::unique_ptr<TH1>(cloned_hist);
}

std::string StackedHistogram::format_yield(double value, int precision) const {
  std::ostringstream ss;
  ss.setf(std::ios::fixed);
  ss.precision(precision);
  ss << value;
  return ss.str();
}

void StackedHistogram::set_global_style() const {
  faint::plot::apply_plot_style();
}

void StackedHistogram::draw_cuts(double max_y, double x_min, double x_max) {
  if (cuts_.empty()) {
    return;
  }

  const double y_for_arrow = max_y * 0.85;
  const double arrow_length = (x_max - x_min) * 0.04;

  for (const auto& cut : cuts_) {
    auto line = std::make_unique<TLine>(cut.threshold, 0.0, cut.threshold, max_y);
    line->SetLineColor(cut.color);
    line->SetLineWidth(2);
    line->SetLineStyle(kDashed);
    line->Draw("same");

    double x_end = cut.direction == CutDirection::kGreaterThan
                       ? cut.threshold + arrow_length
                       : cut.threshold - arrow_length;

    auto arrow = std::make_unique<TArrow>(cut.threshold, y_for_arrow, x_end,
                                          y_for_arrow, 0.02, ">");
    arrow->SetLineColor(cut.color);
    arrow->SetFillColor(cut.color);
    arrow->SetLineWidth(2);
    arrow->Draw("same");

    overlays_.push_back(std::unique_ptr<TObject>(line.release()));
    overlays_.push_back(std::unique_ptr<TObject>(arrow.release()));

    if (!cut.label.empty()) {
      auto text = std::make_unique<TLatex>(cut.threshold, max_y * 1.02,
                                          cut.label.c_str());
      text->SetTextAlign(21);
      text->SetTextFont(42);
      text->SetTextSize(0.04);
      text->Draw();
      overlays_.push_back(std::unique_ptr<TObject>(text.release()));
    }
  }
}

}  // namespace plot
}  // namespace faint
