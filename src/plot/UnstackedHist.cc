#include "rarexsec/plot/UnstackedHist.hh"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TH1D.h>
#include <TLatex.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <utility>
#include <vector>

#include "rarexsec/Plotter.hh"                // H1Spec, Options, Plotter::fmt_commas/sanitise
#include "rarexsec/plot/Channels.hh"          // channel colors/labels/order
#include "rarexsec/proc/Selection.hh"              // selection::apply
#include "rarexsec/Hub.hh"                    // Entry

namespace {
inline void normalise_pdf(TH1D& h) {
  const double area = h.Integral("width");
  if (area > 0.0) h.Scale(1.0 / area);
}
} // namespace

namespace rarexsec::plot {

UnstackedHist::UnstackedHist(H1Spec spec,
                             Options opt,
                             std::vector<const Entry*> mc,
                             std::vector<const Entry*> data,
                             bool normalize_to_pdf,
                             int line_width)
  : spec_(std::move(spec))
  , opt_(std::move(opt))
  , mc_(std::move(mc))
  , data_(std::move(data))
  , normalize_to_pdf_(normalize_to_pdf)
  , line_width_(line_width)
  , plot_name_(rarexsec::plot::Plotter::sanitise(spec_.id))
  , output_directory_(opt_.out_dir) {}

void UnstackedHist::setup_pads(TCanvas& c, TPad*& p_main, TPad*& p_legend) const {
  c.cd();
  p_main = nullptr; p_legend = nullptr;
  const double split = std::clamp(opt_.legend_split, 0.60, 0.95);

  if (opt_.legend_on_top) {
    p_main   = new TPad("pad_main",   "pad_main",   0., 0.00, 1., split);
    p_legend = new TPad("pad_legend", "pad_legend", 0., split, 1., 1.00);

    p_main  ->SetTopMargin(0.01);  p_main  ->SetBottomMargin(0.12);
    p_main  ->SetLeftMargin(0.14); p_main  ->SetRightMargin(0.05);  // better y-axis margin
    p_legend->SetTopMargin(0.05);  p_legend->SetBottomMargin(0.01);
    p_legend->SetLeftMargin(0.02); p_legend->SetRightMargin(0.02);

    if (opt_.use_log_y) p_main->SetLogy();

    p_main->Draw(); p_legend->Draw();
  } else {
    p_main  = new TPad("pad_main","pad_main", 0.,0.,1.,1.);
    p_main ->SetTopMargin(0.06);  p_main ->SetBottomMargin(0.12);
    p_main ->SetLeftMargin(0.14); p_main ->SetRightMargin(0.05);   // better y-axis margin
    if (opt_.use_log_y) p_main->SetLogy();
    p_main->Draw();
  }
}

void UnstackedHist::build_histograms() {
  mc_ch_hists_.clear();
  data_hist_.reset();
  chan_order_.clear();

  // Book per-channel histograms across all MC entries
  std::map<int, std::vector<ROOT::RDF::RResultPtr<TH1D>>> booked_mc;
  const auto& channels = rarexsec::plot::Channels::mc_keys();

  for (size_t ie = 0; ie < mc_.size(); ++ie) {
    const Entry* e = mc_[ie];
    if (!e) continue;

    auto n0 = selection::apply(e->rnode(), spec_.sel, *e);
    auto n  = (spec_.expr.empty() ? n0 : n0.Define("_rx_expr_", spec_.expr));
    const std::string var = spec_.expr.empty() ? spec_.id : "_rx_expr_";

    for (int ch : channels) {
      auto nf = n.Filter([ch](int c){ return c == ch; }, {"analysis_channels"});
      auto h  = nf.Histo1D(spec_.model("_mc_ch"+std::to_string(ch)+"_src"+std::to_string(ie)), var, spec_.weight);
      booked_mc[ch].push_back(h);
    }
  }

  // Sum per-channel across entries; order by pre-normalization yield
  std::vector<std::pair<int,double>> yields;
  std::map<int, std::unique_ptr<TH1D>> sum_by_channel;

  for (int ch : channels) {
    auto it = booked_mc.find(ch);
    if (it == booked_mc.end() || it->second.empty()) continue;

    std::unique_ptr<TH1D> sum;
    for (auto& rr : it->second) {
      const TH1D& h = rr.GetValue();
      if (!sum) {
        sum.reset(static_cast<TH1D*>(h.Clone((spec_.id+"_sum_ch"+std::to_string(ch)).c_str())));
        sum->SetDirectory(nullptr);
      } else {
        sum->Add(&h);
      }
    }
    if (sum) {
      const double y = sum->Integral();
      yields.emplace_back(ch, y);
      sum_by_channel.emplace(ch, std::move(sum));
    }
  }

  std::stable_sort(yields.begin(), yields.end(), [](const auto& a, const auto& b){
    if (a.second == b.second) return a.first < b.first;
    return a.second > b.second;
  });

  // Style: lines by channel color; optional PDF normalization
  for (auto& [ch, /*y*/_] : yields) {
    auto it = sum_by_channel.find(ch);
    if (it == sum_by_channel.end()) continue;
    auto& h = it->second;

    if (normalize_to_pdf_) normalise_pdf(*h);

    h->SetFillStyle(0);
    h->SetLineColor(rarexsec::plot::Channels::color(ch));
    h->SetLineWidth(line_width_);

    mc_ch_hists_.push_back(std::move(h));
    chan_order_.push_back(ch);
  }

  // Optional data overlay (pooled)
  if (!data_.empty()) {
    std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
    for (size_t ie = 0; ie < data_.size(); ++ie) {
      const Entry* e = data_[ie];
      if (!e) continue;
      auto n0 = selection::apply(e->rnode(), spec_.sel, *e);
      auto n  = (spec_.expr.empty() ? n0 : n0.Define("_rx_expr_", spec_.expr));
      const std::string var = spec_.expr.empty() ? spec_.id : "_rx_expr_";
      parts.push_back(n.Histo1D(spec_.model("_data_src"+std::to_string(ie)), var));
    }
    for (auto& rr : parts) {
      const TH1D& h = rr.GetValue();
      if (!data_hist_) {
        data_hist_.reset(static_cast<TH1D*>(h.Clone((spec_.id+"_data").c_str())));
        data_hist_->SetDirectory(nullptr);
      } else {
        data_hist_->Add(&h);
      }
    }
    if (data_hist_) {
      if (normalize_to_pdf_) normalise_pdf(*data_hist_);
      data_hist_->SetMarkerStyle(kFullCircle);
      data_hist_->SetMarkerSize(0.9);
      data_hist_->SetLineColor(kBlack);
      data_hist_->SetFillStyle(0);
    }
  }
}

void UnstackedHist::draw_curves(TPad* p_main, double& max_y) {
  if (!p_main) return;
  p_main->cd();

  // find max
  max_y = 0.;
  for (auto& h : mc_ch_hists_) max_y = std::max(max_y, h->GetMaximum());
  if (data_hist_) max_y = std::max(max_y, data_hist_->GetMaximum());
  if (opt_.y_max > 0) max_y = opt_.y_max;

  if (mc_ch_hists_.empty()) return;
  TH1D* frame = mc_ch_hists_.front().get();

  if (spec_.xmin < spec_.xmax) frame->GetXaxis()->SetLimits(spec_.xmin, spec_.xmax);

  // Title: only axis titles, no main title.
  frame->SetTitle((std::string(";") + opt_.x_title + ";" + (normalize_to_pdf_ ? "Probability density" : opt_.y_title)).c_str());

  // NEW: readable numeric scale (no ×10^n)
  frame->GetXaxis()->SetNoExponent(true);
  frame->GetYaxis()->SetNoExponent(true);
  frame->GetXaxis()->SetMaxDigits(4);
  frame->GetYaxis()->SetMaxDigits(4);

  frame->SetMaximum(max_y * (opt_.use_log_y ? 10. : 1.3));
  frame->SetMinimum(opt_.use_log_y ? 0.1 : opt_.y_min);

  frame->Draw("HIST");
  for (size_t i = 1; i < mc_ch_hists_.size(); ++i) mc_ch_hists_[i]->Draw("HIST SAME");
  if (data_hist_) data_hist_->Draw("E1 SAME");
}

void UnstackedHist::draw_legend(TPad* p) {
  if (!p) return;
  p->cd();
  legend_ = std::make_unique<TLegend>(0.12, 0.0, 0.95, 0.75);
  auto* leg = legend_.get();
  if (opt_.legend_on_top) {
    // use (nearly) the whole top legend pad and multiple columns like StackedHist
    leg->SetX1NDC(0.02); leg->SetY1NDC(0.05);
    leg->SetX2NDC(0.98); leg->SetY2NDC(0.95);
  } else {
    leg->SetX1NDC(opt_.leg_x1); leg->SetY1NDC(opt_.leg_y1);
    leg->SetX2NDC(opt_.leg_x2); leg->SetY2NDC(opt_.leg_y2);
  }
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);

  // multi-column layout when there are many entries
  int n_entries = static_cast<int>(mc_ch_hists_.size()) + (data_hist_ ? 1 : 0);
  if (n_entries > 0) leg->SetNColumns(n_entries > 4 ? 3 : 2);

  legend_proxies_.clear();
  for (size_t i = 0; i < mc_ch_hists_.size(); ++i) {
    int ch = chan_order_.at(i);
    auto proxy = std::unique_ptr<TH1D>(static_cast<TH1D*>(
      mc_ch_hists_[i]->Clone((spec_.id+"_leg_ch"+std::to_string(ch)).c_str())));
    proxy->SetDirectory(nullptr);
    proxy->Reset("ICES");

    proxy->SetFillStyle(0);
    proxy->SetLineColor(mc_ch_hists_[i]->GetLineColor());
    proxy->SetLineWidth(line_width_);

    std::string label = rarexsec::plot::Channels::label(ch);
    if (opt_.annotate_numbers) {
      const double sum = mc_ch_hists_[i]->Integral();
      label += " : " + rarexsec::plot::Plotter::fmt_commas(sum, 2);
    }
    leg->AddEntry(proxy.get(), label.c_str(), "l");
    legend_proxies_.push_back(std::move(proxy));
  }
  if (data_hist_) leg->AddEntry(data_hist_.get(), "Data", "lep");
  leg->Draw();
}

void UnstackedHist::draw_watermark(TPad* p_main) const {
  if (!p_main) return;
  p_main->cd();
  TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.04);
  tl.DrawLatex(0.14, 0.92,
     (std::string("#muBooNE Simulation – ")
      + (opt_.analysis_region_label.empty() ? "Empty Selection" : opt_.analysis_region_label)).c_str());
}

void UnstackedHist::draw(TCanvas& canvas) {
  build_histograms();
  TPad *p_main = nullptr, *p_legend = nullptr;
  setup_pads(canvas, p_main, p_legend);
  double max_y = 1.;
  draw_curves(p_main, max_y);
  draw_watermark(p_main);
  draw_legend(p_legend ? p_legend : p_main);
  if (p_main) p_main->RedrawAxis();
  canvas.Update();
}

void UnstackedHist::draw_and_save(const std::string& image_format) {
  std::filesystem::create_directories(output_directory_);
  TCanvas canvas(plot_name_.c_str(), plot_name_.c_str(), 900, 700);
  draw(canvas);
  const std::string fmt = image_format.empty() ? "png" : image_format;
  canvas.SaveAs((output_directory_ + "/" + plot_name_ + "." + fmt).c_str());
}

} // namespace rarexsec::plot
