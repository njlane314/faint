#include "rarexsec/plot/StackedHist.hh"
#include "ROOT/RDFHelpers.hxx"
#include "TArrow.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMatrixDSym.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include "rarexsec/plot/Channels.hh"

namespace {

static void apply_total_errors(TH1D& h, const TMatrixDSym* cov, const std::vector<double>* syst_bin) {
    const int nb = h.GetNbinsX();
    for (int i = 1; i <= nb; ++i) {
        const double stat = h.GetBinError(i);
        double syst = 0.0;
        if (cov && i - 1 < cov->GetNrows()) {
            syst = std::sqrt((*cov)(i - 1, i - 1));
        } else if (syst_bin && i - 1 < static_cast<int>(syst_bin->size())) {
            syst = std::max(0.0, (*syst_bin)[i - 1]);
        }
        const double tot = std::sqrt(stat * stat + syst * syst);
        h.SetBinError(i, tot);
    }
}

}

rarexsec::plot::StackedHist::StackedHist(H1Spec spec,
                              Options opt,
                              std::vector<const Entry*> mc,
                              std::vector<const Entry*> data)
: spec_(std::move(spec))
, opt_(std::move(opt))
, mc_(std::move(mc))
, data_(std::move(data))
, plot_name_(rarexsec::plot::Plotter::sanitise(spec_.id))
, output_directory_(opt_.out_dir) {}

void rarexsec::plot::StackedHist::setup_pads(TCanvas& c, TPad*& p_main, TPad*& p_ratio, TPad*& p_legend) const {
    c.cd();
    p_main = nullptr;
    p_ratio = nullptr;
    p_legend = nullptr;

    const double split = std::clamp(opt_.legend_split, 0.60, 0.95);

    if (opt_.legend_on_top) {
        if (want_ratio()) {
            p_ratio  = new TPad("pad_ratio",  "pad_ratio",   0., 0.00, 1., 0.30);
            p_main   = new TPad("pad_main",   "pad_main",    0., 0.30, 1., split);
            p_legend = new TPad("pad_legend", "pad_legend",  0., split, 1., 1.00);

            p_main ->SetTopMargin(0.02);   p_main ->SetBottomMargin(0.02);
            p_main ->SetLeftMargin(0.12);  p_main ->SetRightMargin(0.05);

            p_ratio->SetTopMargin(0.05);   p_ratio->SetBottomMargin(0.35);
            p_ratio->SetLeftMargin(0.12);  p_ratio->SetRightMargin(0.05);

            p_legend->SetTopMargin(0.05);  p_legend->SetBottomMargin(0.01);
            p_legend->SetLeftMargin(0.02); p_legend->SetRightMargin(0.02);
        } else {
            p_main   = new TPad("pad_main",   "pad_main",   0., 0.00, 1., split);
            p_legend = new TPad("pad_legend", "pad_legend", 0., split, 1., 1.00);

            p_main ->SetTopMargin(0.01);  p_main ->SetBottomMargin(0.12);
            p_main ->SetLeftMargin(0.12); p_main ->SetRightMargin(0.05);

            p_legend->SetTopMargin(0.05); p_legend->SetBottomMargin(0.01);
            p_legend->SetLeftMargin(0.02); p_legend->SetRightMargin(0.02);
        }
        if (opt_.use_log_y && p_main) p_main->SetLogy();
        if (p_ratio)  p_ratio->Draw();
        if (p_main)   p_main->Draw();
        if (p_legend) p_legend->Draw();
    } else {
        if (want_ratio()) {
            p_main  = new TPad("pad_main","pad_main", 0.,0.30,1.,1.);
            p_ratio = new TPad("pad_ratio","pad_ratio",0.,0.,  1.,0.30);
            p_main ->SetTopMargin(0.06);  p_main ->SetBottomMargin(0.02);
            p_main ->SetLeftMargin(0.12); p_main ->SetRightMargin(0.05);
            p_ratio->SetTopMargin(0.05);  p_ratio->SetBottomMargin(0.35);
            p_ratio->SetLeftMargin(0.12); p_ratio->SetRightMargin(0.05);
            if (opt_.use_log_y) p_main->SetLogy();
            p_ratio->Draw(); p_main->Draw();
        } else {
            p_main  = new TPad("pad_main","pad_main", 0.,0.,1.,1.);
            p_main ->SetTopMargin(0.06);  p_main ->SetBottomMargin(0.12);
            p_main ->SetLeftMargin(0.12); p_main ->SetRightMargin(0.05);
            if (opt_.use_log_y) p_main->SetLogy();
            p_main->Draw();
        }
    }
}

void rarexsec::plot::StackedHist::build_histograms() {
    const auto axes = spec_.axis_title();
    stack_ = std::make_unique<THStack>((spec_.id + "_stack").c_str(), axes.c_str());
    mc_ch_hists_.clear();
    mc_total_.reset();
    data_hist_.reset();
    sig_hist_.reset();
    signal_scale_ = 1.0;
    std::map<int, std::vector<ROOT::RDF::RResultPtr<TH1D>>> booked;
    const auto& channels = rarexsec::plot::Channels::mc_keys();

    auto rebin_hist = [&](std::unique_ptr<TH1D>& hist, const std::string& suffix) {
        if (!hist || opt_.rebin_edges.size() < 2) return;
        const int nb = static_cast<int>(opt_.rebin_edges.size()) - 1;
        std::string newname = hist->GetName();
        newname += suffix;
        auto* rebinned = hist->Rebin(nb, newname.c_str(), opt_.rebin_edges.data());
        hist.reset(static_cast<TH1D*>(rebinned));
        hist->SetDirectory(nullptr);
    };

    for (size_t ie = 0; ie < mc_.size(); ++ie) {
        const Entry* e = mc_[ie];
        if (!e) continue;
        auto n0 = selection::apply(e->rnode(), spec_.sel, *e);
        auto n  = (spec_.expr.empty() ? n0 : n0.Define("_rx_expr_", spec_.expr));
        const std::string var = spec_.expr.empty() ? spec_.id : "_rx_expr_";
        for (int ch : channels) {
            auto nf = n.Filter([ch](int c){ return c==ch; }, {"analysis_channels"});
            auto h  = nf.Histo1D(spec_.model("_mc_ch"+std::to_string(ch)+"_src"+std::to_string(ie)), var, spec_.weight);
            booked[ch].push_back(h);
        }
    }

    std::vector<int> order;
    std::map<int, std::unique_ptr<TH1D>> sum_by_channel;
    std::vector<std::pair<int,double>> yields;

    for (int ch : channels) {
        auto it = booked.find(ch);
        if (it == booked.end() || it->second.empty()) continue;
        std::unique_ptr<TH1D> sum;
        for (auto& rr : it->second) {
            const TH1D& h = rr.GetValue();
            if (!sum) {
                sum.reset(static_cast<TH1D*>(h.Clone((spec_.id+"_mc_sum_ch"+std::to_string(ch)).c_str())));
                sum->SetDirectory(nullptr);
            } else {
                sum->Add(&h);
            }
        }
        if (sum) {
            rebin_hist(sum, "_rebin");
            double y = sum->Integral();
            yields.emplace_back(ch, y);
            sum_by_channel.emplace(ch, std::move(sum));
        }
    }

    std::stable_sort(yields.begin(), yields.end(), [](const auto& a, const auto& b){
        if (a.second == b.second) return a.first < b.first;
        return a.second > b.second;
    });

    chan_order_.clear();
    for (auto& cy : yields) {
        int ch = cy.first;
        auto it = sum_by_channel.find(ch);
        if (it == sum_by_channel.end()) continue;
        auto& sum = it->second;
        sum->SetFillColor(rarexsec::plot::Channels::color(ch));
        sum->SetFillStyle(rarexsec::plot::Channels::fill_style(ch));
        sum->SetLineColor(kBlack);
        sum->SetLineWidth(1);
        stack_->Add(sum.get(), "HIST");
        mc_ch_hists_.push_back(std::move(sum));
        chan_order_.push_back(ch);
    }

    for (auto& uptr : mc_ch_hists_) {
        if (!mc_total_) {
            mc_total_.reset(static_cast<TH1D*>(uptr->Clone((spec_.id+"_mc_total").c_str())));
            mc_total_->SetDirectory(nullptr);
        } else {
            mc_total_->Add(uptr.get());
        }
    }

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
            rebin_hist(data_hist_, "_rebin");
            data_hist_->SetMarkerStyle(kFullCircle);
            data_hist_->SetMarkerSize(0.9);
            data_hist_->SetLineColor(kBlack);
            data_hist_->SetFillStyle(0);
        }
    }

    if (opt_.overlay_signal && !opt_.signal_channels.empty() && !mc_ch_hists_.empty()) {
        double tot_sum = mc_total_ ? mc_total_->Integral() : 0.0;
        auto sig = std::make_unique<TH1D>(*mc_ch_hists_.front());
        sig->Reset();
        for (size_t i = 0; i < mc_ch_hists_.size(); ++i) {
            int ch = chan_order_.at(i);
            if (std::find(opt_.signal_channels.begin(), opt_.signal_channels.end(), ch) != opt_.signal_channels.end()) {
                sig->Add(mc_ch_hists_[i].get());
            }
        }
        double sig_sum = sig->Integral();
        if (sig_sum > 0.0 && tot_sum > 0.0) {
            signal_scale_ = tot_sum / sig_sum;
            sig->Scale(signal_scale_);
        }
        sig->SetLineColor(kGreen+2);
        sig->SetLineStyle(kDashed);
        sig->SetLineWidth(2);
        sig->SetFillStyle(0);
        sig_hist_ = std::move(sig);
    }
}

void rarexsec::plot::StackedHist::draw_stack_and_unc(TPad* p_main, double& max_y) {
    if (!p_main) return;
    p_main->cd();
    stack_->Draw("HIST");
    TH1* frame = stack_->GetHistogram();
    if (frame) {
        frame->SetLineWidth(2);
    }
    if (frame && spec_.xmin < spec_.xmax) frame->GetXaxis()->SetLimits(spec_.xmin, spec_.xmax);
    if (frame) {
        frame->GetXaxis()->SetNdivisions(510);
        frame->GetXaxis()->SetTickLength(0.02);
        if (!opt_.x_title.empty()) frame->GetXaxis()->SetTitle(opt_.x_title.c_str());
        if (!opt_.y_title.empty()) frame->GetYaxis()->SetTitle(opt_.y_title.c_str());
    }
    if (mc_total_) {
        if (opt_.total_cov || !opt_.syst_bin.empty()) {
            apply_total_errors(*mc_total_, opt_.total_cov.get(),
                               opt_.syst_bin.empty() ? nullptr : &opt_.syst_bin);
        }

        max_y = mc_total_->GetMaximum() + mc_total_->GetBinError(mc_total_->GetMaximumBin());
        if (opt_.y_max > 0) max_y = opt_.y_max;

        stack_->SetMaximum(max_y * (opt_.use_log_y ? 10. : 1.3));
        stack_->SetMinimum(opt_.use_log_y ? 0.1 : opt_.y_min);

        auto* h = static_cast<TH1D*>(mc_total_->Clone((spec_.id+"_mc_totband").c_str()));
        h->SetDirectory(nullptr);
        h->SetFillColor(kBlack);
        h->SetFillStyle(3004);
        h->SetMarkerSize(0);
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
        h->Draw("E2 SAME");
        h->Draw("E1 SAME");
    }
    if (sig_hist_)  sig_hist_->Draw("HIST SAME");
    if (data_hist_) data_hist_->Draw("E1 SAME");
}

void rarexsec::plot::StackedHist::draw_ratio(TPad* p_ratio) {
    if (!p_ratio || !data_hist_ || !mc_total_) return;
    p_ratio->cd();

    auto ratio = std::unique_ptr<TH1D>(static_cast<TH1D*>(
        data_hist_->Clone((spec_.id+"_ratio").c_str())));
    ratio->SetDirectory(nullptr);
    ratio->Divide(mc_total_.get());
    ratio->SetTitle("; ;Data / MC");
    ratio->SetMaximum(1.99);
    ratio->SetMinimum(0.01);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetXaxis()->SetLabelSize(0.10);
    ratio->GetYaxis()->SetLabelSize(0.10);
    ratio->GetYaxis()->SetTitleSize(0.10);
    ratio->GetYaxis()->SetTitleOffset(0.4);
    ratio->GetYaxis()->SetTitle("Data / MC");
    ratio->GetXaxis()->SetTitle(opt_.x_title.empty() ? data_hist_->GetXaxis()->GetTitle() : opt_.x_title.c_str());

    ratio->Draw("E1");

    std::unique_ptr<TH1D> band;
    if (opt_.show_ratio_band) {
        band.reset(static_cast<TH1D*>(mc_total_->Clone((spec_.id+"_ratio_band").c_str())));
        band->SetDirectory(nullptr);
        const int nb = band->GetNbinsX();
        for (int i = 1; i <= nb; ++i) {
            const double m  = mc_total_->GetBinContent(i);
            const double em = mc_total_->GetBinError(i);
            band->SetBinContent(i, 1.0);
            band->SetBinError(i, (m > 0 ? em / m : 0.0));
        }
        band->SetFillColor(kBlack);
        band->SetFillStyle(3004);
        band->SetLineColor(kBlack);
        band->SetMarkerSize(0);
        band->Draw("E2 SAME");
    }

    ratio->Draw("E1 SAME");
}

void rarexsec::plot::StackedHist::draw_legend(TPad* p) {
    if (!p) return;
    p->cd();
    legend_ = std::make_unique<TLegend>(0.12, 0.0, 0.95, 0.75);
    auto* leg = legend_.get();
    if (!opt_.legend_on_top) {
        leg->SetX1NDC(opt_.leg_x1); leg->SetY1NDC(opt_.leg_y1);
        leg->SetX2NDC(opt_.leg_x2); leg->SetY2NDC(opt_.leg_y2);
    }
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);

    int n_entries = static_cast<int>(mc_ch_hists_.size());
    if (mc_total_) ++n_entries;
    if (sig_hist_) ++n_entries;
    if (data_hist_) ++n_entries;
    if (n_entries > 0) leg->SetNColumns(n_entries > 4 ? 3 : 2);

    legend_proxies_.clear();

    for (size_t i = 0; i < mc_ch_hists_.size(); ++i) {
        int ch = chan_order_.at(i);
        double sum = mc_ch_hists_[i]->Integral();
        std::string label = rarexsec::plot::Channels::label(ch);
        if (opt_.annotate_numbers) {
            label += " : " + rarexsec::plot::Plotter::fmt_commas(sum, 2);
        }
        auto proxy = std::unique_ptr<TH1D>(static_cast<TH1D*>(
            mc_ch_hists_[i]->Clone((spec_.id+"_leg_ch"+std::to_string(ch)).c_str())));
        proxy->SetDirectory(nullptr);
        proxy->Reset("ICES");

        auto* entry = leg->AddEntry(proxy.get(), label.c_str(), "f");
        leg->SetEntrySeparation(0.01);
        legend_proxies_.push_back(std::move(proxy));
        (void)entry;
    }

    if (mc_total_) {
        auto proxy = std::unique_ptr<TH1D>(static_cast<TH1D*>(
            mc_total_->Clone((spec_.id+"_leg_unc").c_str())));
        proxy->SetDirectory(nullptr);
        proxy->Reset("ICES");
        proxy->SetFillColor(kBlack);
        proxy->SetFillStyle(3004);
        proxy->SetLineColor(kBlack);
        proxy->SetLineWidth(1);
        leg->AddEntry(proxy.get(), "Stat. #oplus Syst. Unc.", "f");
        legend_proxies_.push_back(std::move(proxy));
    }

    if (sig_hist_) {
        std::string sig_label = "Signal";
        if (signal_scale_ != 1.0) {
            sig_label += " (x" + rarexsec::plot::Plotter::fmt_commas(signal_scale_, 2) + ")";
        }
        leg->AddEntry(sig_hist_.get(), sig_label.c_str(), "l");
    }

    if (data_hist_) {
        leg->AddEntry(data_hist_.get(), "Data", "lep");
    }

    leg->Draw();
}

void rarexsec::plot::StackedHist::draw_cuts(TPad* p, double max_y) {
    if (!opt_.show_cuts || opt_.cuts.empty()) return;
    p->cd();
    TH1* frame = stack_->GetHistogram();
    if (!frame) return;
    const double y = max_y * 0.85;
    const double xmin = frame->GetXaxis()->GetXmin();
    const double xmax = frame->GetXaxis()->GetXmax();
    const double xr = xmax - xmin;
    const double alen = xr * 0.04;
    for (const auto& c : opt_.cuts) {
        auto* line = new TLine(c.x, 0., c.x, max_y * 1.3);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->SetLineStyle(kDashed);
        line->Draw("same");
        const double x1 = c.x;
        const double x2 = (c.dir == CutDir::GreaterThan) ? c.x + alen : c.x - alen;
        auto* arr = new TArrow(x1, y, x2, y, 0.025, ">");
        arr->SetLineColor(kRed);
        arr->SetFillColor(kRed);
        arr->SetLineWidth(2);
        arr->Draw("same");
    }
}

void rarexsec::plot::StackedHist::draw_watermark(TPad* p, double total_mc) const {
    if (!p) return;
    p->cd();
    TLatex lt;
    lt.SetNDC();
    const double x = p->GetLeftMargin() + 0.03;
    double y = 1 - p->GetTopMargin() - 0.03;

    const std::string title = opt_.watermark_title.empty()
                                ? std::string("rarexsec, Preliminary")
                                : opt_.watermark_title;

    lt.SetTextAlign(13);
    lt.SetTextFont(62);
    lt.SetTextSize(0.05);
    lt.DrawLatex(x, y, ("#bf{" + title + "}").c_str());

    std::vector<std::string> lines = opt_.watermark_lines;
    if (lines.empty()) {
        std::string bl = opt_.beamline.empty() ? "N/A" : opt_.beamline;
        std::string runs = opt_.periods.empty() ? "N/A" : [&]{
            std::string s;
            for (size_t i = 0; i < opt_.periods.size(); ++i) {
                s += opt_.periods[i];
                if (i + 1 < opt_.periods.size()) s += ", ";
            }
            return s;
        }();
        lines.push_back("Beamline, Periods: " + bl + ", " + runs);
        lines.push_back("Total MC: " + rarexsec::plot::Plotter::fmt_commas(total_mc, 2) + " events");
    }

    lt.SetTextFont(42);
    const double line_size = 0.05 * 0.8;
    lt.SetTextSize(line_size);
    const double step = line_size * 1.2;
    for (const auto& line : lines) {
        y -= step;
        if (y < 0) break;
        lt.DrawLatex(x, y, line.c_str());
    }
}

void rarexsec::plot::StackedHist::draw(TCanvas& canvas) {
    build_histograms();
    TPad *p_main = nullptr, *p_ratio = nullptr, *p_legend = nullptr;
    setup_pads(canvas, p_main, p_ratio, p_legend);
    double max_y = 1.;
    draw_stack_and_unc(p_main, max_y);
    draw_cuts(p_main, max_y);
    draw_watermark(p_main, mc_total_ ? mc_total_->Integral() : 0.0);
    draw_legend(p_legend ? p_legend : p_main);
    if (want_ratio()) draw_ratio(p_ratio);
    if (p_main) p_main->RedrawAxis();
    canvas.Update();
}

void rarexsec::plot::StackedHist::draw_and_save(const std::string& image_format) {
    std::filesystem::create_directories(output_directory_);
    TCanvas canvas(plot_name_.c_str(), plot_name_.c_str(), 800, 600);
    draw(canvas);
    const std::string fmt = image_format.empty() ? "png" : image_format;
    canvas.SaveAs((output_directory_ + "/" + plot_name_ + "." + fmt).c_str());
}