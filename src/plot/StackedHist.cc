#include "rarexsec/plot/StackedHist.hh"
#include "ROOT/RDFHelpers.hxx"
#include "TArrow.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include <algorithm>
#include <filesystem>
#include <map>
#include "rarexsec/plot/Channels.hh"

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

void rarexsec::plot::StackedHist::setup_pads(TCanvas& c, TPad*& p_main, TPad*& p_ratio) const {
    c.cd();
    if (want_ratio()) {
        p_main  = new TPad("pad_main","pad_main", 0.,0.30,1.,1.);
        p_ratio = new TPad("pad_ratio","pad_ratio",0.,0.,  1.,0.30);
        p_main ->SetTopMargin(0.06);  p_main ->SetBottomMargin(0.02);
        p_main ->SetLeftMargin(0.12); p_main ->SetRightMargin(0.05);
        p_ratio->SetTopMargin(0.05);  p_ratio->SetBottomMargin(0.35);
        p_ratio->SetLeftMargin(0.12); p_ratio->SetRightMargin(0.05);
        if (opt_.use_log_y) p_main->SetLogy();
        p_main->Draw(); p_ratio->Draw();
    } else {
        p_main  = new TPad("pad_main","pad_main", 0.,0.,1.,1.);
        p_ratio = nullptr;
        p_main ->SetTopMargin(0.06);  p_main ->SetBottomMargin(0.12);
        p_main ->SetLeftMargin(0.12); p_main ->SetRightMargin(0.05);
        if (opt_.use_log_y) p_main->SetLogy();
        p_main->Draw();
    }
}

void rarexsec::plot::StackedHist::build_histograms() {
    stack_ = std::make_unique<THStack>((spec_.id + "_stack").c_str(), spec_.title.c_str());
    std::map<int, std::vector<ROOT::RDF::RResultPtr<TH1D>>> booked;
    const auto& channels = rarexsec::plot::Channels::mc_keys();

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
            data_hist_->SetMarkerStyle(kFullCircle);
            data_hist_->SetMarkerSize(0.9);
            data_hist_->SetLineColor(kBlack);
            data_hist_->SetFillStyle(0);
        }
    }

    if (opt_.overlay_signal && !mc_ch_hists_.empty()) {
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
        if (sig_sum > 0.0 && tot_sum > 0.0) sig->Scale(tot_sum / sig_sum);
        sig->SetLineColor(kGreen+2);
        sig->SetLineStyle(kDashed);
        sig->SetLineWidth(2);
        sig->SetFillStyle(0);
        sig_hist_ = std::move(sig);
    }
}

void rarexsec::plot::StackedHist::draw_stack_and_unc(TPad* p_main, double& max_y) {
    p_main->cd();
    stack_->Draw("HIST");
    TH1* frame = stack_->GetHistogram();
    if (frame && spec_.xmin < spec_.xmax) frame->GetXaxis()->SetLimits(spec_.xmin, spec_.xmax);
    if (frame) {
        frame->GetXaxis()->SetNdivisions(510);
        frame->GetXaxis()->SetTickLength(0.02);
    }
    if (mc_total_) {
        max_y = mc_total_->GetMaximum() + mc_total_->GetBinError(mc_total_->GetMaximumBin());
        if (opt_.y_max > 0) max_y = opt_.y_max;
        stack_->SetMaximum(max_y * (opt_.use_log_y ? 10. : 1.3));
        stack_->SetMinimum(opt_.use_log_y ? 0.1 : opt_.y_min);
        auto* h = static_cast<TH1D*>(mc_total_->Clone((spec_.id+"_mc_statband").c_str()));
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
    auto ratio = std::unique_ptr<TH1D>(static_cast<TH1D*>(data_hist_->Clone((spec_.id+"_ratio").c_str())));
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
    ratio->Draw("E1");
}

void rarexsec::plot::StackedHist::draw_legend(TPad* p) {
    p->cd();
    auto* leg = new TLegend(opt_.leg_x1, opt_.leg_y1, opt_.leg_x2, opt_.leg_y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    for (size_t i = 0; i < mc_ch_hists_.size(); ++i) {
        int ch = chan_order_.at(i);
        auto* h = mc_ch_hists_[i].get();
        double sum = h->Integral();
        auto lab = rarexsec::plot::Channels::label(ch);
        if (opt_.annotate_numbers) leg->AddEntry(h, (lab + " : " + rarexsec::plot::Plotter::fmt_commas(sum, 2)).c_str(), "f");
        else                       leg->AddEntry(h, lab.c_str(), "f");
    }
    if (mc_total_) {
        auto* h_unc = static_cast<TH1D*>(mc_total_->Clone());
        h_unc->SetFillColor(kBlack);
        h_unc->SetFillStyle(3004);
        h_unc->SetLineColor(kBlack);
        h_unc->SetLineWidth(1);
        h_unc->SetMarkerSize(0);
        leg->AddEntry(h_unc, "Stat. Unc.", "f");
    }
    if (sig_hist_)  leg->AddEntry(sig_hist_.get(), "Signal (scaled)", "l");
    if (data_hist_) leg->AddEntry(data_hist_.get(), "Data", "lep");
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
    p->cd();
    const std::string bl = opt_.beamline.empty() ? "N/A" : opt_.beamline;
    std::string runs;
    if (!opt_.periods.empty()) {
        for (size_t i = 0; i < opt_.periods.size(); ++i) {
            runs += opt_.periods[i];
            if (i + 1 < opt_.periods.size()) runs += ", ";
        }
    } else {
        runs = "N/A";
    }
    TLatex lt;
    lt.SetNDC();
    lt.SetTextAlign(33);
    lt.SetTextFont(62);
    lt.SetTextSize(0.05);
    lt.DrawLatex(1 - p->GetRightMargin() - 0.03, 1 - p->GetTopMargin() - 0.03, "#bf{rarexsec, Preliminary}");
    lt.SetTextFont(42);
    lt.SetTextSize(0.05 * 0.8);
    lt.DrawLatex(1 - p->GetRightMargin() - 0.03, 1 - p->GetTopMargin() - 0.09, (std::string("Beamline, Periods: ") + bl + ", " + runs).c_str());
    lt.DrawLatex(1 - p->GetRightMargin() - 0.03, 1 - p->GetTopMargin() - 0.15, (std::string("Total MC: ") + rarexsec::plot::Plotter::fmt_commas(total_mc, 2) + " events").c_str());
}

void rarexsec::plot::StackedHist::draw(TCanvas& canvas) {
    build_histograms();
    TPad *p_main = nullptr, *p_ratio = nullptr;
    setup_pads(canvas, p_main, p_ratio);
    double max_y = 1.;
    draw_stack_and_unc(p_main, max_y);
    draw_cuts(p_main, max_y);
    draw_watermark(p_main, mc_total_ ? mc_total_->Integral() : 0.0);
    draw_legend(p_main);
    if (want_ratio()) draw_ratio(p_ratio);
    p_main->RedrawAxis();
    canvas.Update();
}

void rarexsec::plot::StackedHist::draw_and_save(const std::string& image_format) {
    std::filesystem::create_directories(output_directory_);
    TCanvas canvas(plot_name_.c_str(), plot_name_.c_str(), 800, 600);
    draw(canvas);
    const std::string fmt = image_format.empty() ? "png" : image_format;
    canvas.SaveAs((output_directory_ + "/" + plot_name_ + "." + fmt).c_str());
}