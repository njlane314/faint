#include "rarexsec/plot/StackedHist.hh"
#include "ROOT/RDFHelpers.hxx"
#include "TArrow.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TLatex.h"
#include "TLine.h"
#include "TList.h"
#include "TMatrixDSym.h"
#include "rarexsec/plot/Plotter.hh"
#include "rarexsec/plot/Channels.hh"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <map>
#include <sstream>
#include <utility>

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

static std::string selection_label(rarexsec::selection::Preset preset) {
    using rarexsec::selection::Preset;
    switch (preset) {
    case Preset::Trigger:
        return "Trigger Selection";
    case Preset::Slice:
        return "Slice Selection";
    case Preset::Fiducial:
        return "Fiducial Selection";
    case Preset::Topology:
        return "Topology Selection";
    case Preset::Muon:
        return "Muon Selection";
    case Preset::InclusiveMuCC:
        return "Inclusive Muon CC Selection";
    case Preset::Empty:
    default:
        return "Empty Selection";
    }
}

rarexsec::plot::StackedHist::StackedHist(Histogram1DSpec spec,
                                         Options opt,
                                         std::vector<const Entry*> mc,
                                         std::vector<const Entry*> data)
    : spec_(std::move(spec)), opt_(std::move(opt)), mc_(std::move(mc)), data_(std::move(data)), plot_name_(rarexsec::plot::Plotter::sanitise(spec_.id)), output_directory_(opt_.out_dir) {}

void rarexsec::plot::StackedHist::setup_pads(TCanvas& c, TPad*& p_main, TPad*& p_ratio, TPad*& p_legend) const {
    c.cd();
    p_main = nullptr;
    p_ratio = nullptr;
    p_legend = nullptr;

    // Match the reference canvas layout: legend occupies the top 15%
    const double split = 0.85;

    if (opt_.legend_on_top) {
        if (want_ratio()) {
            p_ratio = new TPad("pad_ratio", "pad_ratio", 0., 0.00, 1., 0.30);
            p_main = new TPad("pad_main", "pad_main", 0., 0.30, 1., split);
            p_legend = new TPad("pad_legend", "pad_legend", 0., split, 1., 1.00);

            p_main->SetTopMargin(0.02);
            p_main->SetBottomMargin(0.02);
            p_main->SetLeftMargin(0.12);
            p_main->SetRightMargin(0.05);

            p_ratio->SetTopMargin(0.05);
            p_ratio->SetBottomMargin(0.35);
            p_ratio->SetLeftMargin(0.12);
            p_ratio->SetRightMargin(0.05);

            p_legend->SetTopMargin(0.05);
            p_legend->SetBottomMargin(0.01);
            // Keep full width for the legend, matching the reference implementation
            p_legend->SetLeftMargin(0.00);
            p_legend->SetRightMargin(0.00);
        } else {
            p_main = new TPad("pad_main", "pad_main", 0., 0.00, 1., split);
            p_legend = new TPad("pad_legend", "pad_legend", 0., split, 1., 1.00);

            p_main->SetTopMargin(0.01);
            p_main->SetBottomMargin(0.12);
            p_main->SetLeftMargin(0.12);
            p_main->SetRightMargin(0.05);

            p_legend->SetTopMargin(0.05);
            p_legend->SetBottomMargin(0.01);
            // Keep full width for the legend, matching the reference implementation
            p_legend->SetLeftMargin(0.00);
            p_legend->SetRightMargin(0.00);
        }
        if (opt_.use_log_y && p_main)
            p_main->SetLogy();
        if (p_ratio)
            p_ratio->Draw();
        if (p_main)
            p_main->Draw();
        if (p_legend)
            p_legend->Draw();
    } else {
        if (want_ratio()) {
            p_main = new TPad("pad_main", "pad_main", 0., 0.30, 1., 1.);
            p_ratio = new TPad("pad_ratio", "pad_ratio", 0., 0., 1., 0.30);
            p_main->SetTopMargin(0.06);
            p_main->SetBottomMargin(0.02);
            p_main->SetLeftMargin(0.12);
            p_main->SetRightMargin(0.05);
            p_ratio->SetTopMargin(0.05);
            p_ratio->SetBottomMargin(0.35);
            p_ratio->SetLeftMargin(0.12);
            p_ratio->SetRightMargin(0.05);
            if (opt_.use_log_y)
                p_main->SetLogy();
            p_ratio->Draw();
            p_main->Draw();
        } else {
            p_main = new TPad("pad_main", "pad_main", 0., 0., 1., 1.);
            p_main->SetTopMargin(0.06);
            p_main->SetBottomMargin(0.12);
            p_main->SetLeftMargin(0.12);
            p_main->SetRightMargin(0.05);
            if (opt_.use_log_y)
                p_main->SetLogy();
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

    for (size_t ie = 0; ie < mc_.size(); ++ie) {
        const Entry* e = mc_[ie];
        if (!e)
            continue;
        auto n0 = selection::apply(e->rnode(), spec_.sel, *e);
        auto n = (spec_.expr.empty() ? n0 : n0.Define("_rx_expr_", spec_.expr));
        const std::string var = spec_.expr.empty() ? spec_.id : "_rx_expr_";
        for (int ch : channels) {
            auto nf = n.Filter([ch](int c) { return c == ch; }, {"analysis_channels"});
            auto h = nf.Histo1D(spec_.model("_mc_ch" + std::to_string(ch) + "_src" + std::to_string(ie)), var, spec_.weight);
            booked[ch].push_back(h);
        }
    }

    std::vector<int> order;
    std::map<int, std::unique_ptr<TH1D>> sum_by_channel;
    std::vector<std::pair<int, double>> yields;

    for (int ch : channels) {
        auto it = booked.find(ch);
        if (it == booked.end() || it->second.empty())
            continue;
        std::unique_ptr<TH1D> sum;
        for (auto& rr : it->second) {
            const TH1D& h = rr.GetValue();
            if (!sum) {
                sum.reset(static_cast<TH1D*>(h.Clone((spec_.id + "_mc_sum_ch" + std::to_string(ch)).c_str())));
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

    std::stable_sort(yields.begin(), yields.end(), [](const auto& a, const auto& b) {
        if (a.second == b.second)
            return a.first < b.first;
        return a.second > b.second;
    });

    chan_order_.clear();
    for (auto& cy : yields) {
        int ch = cy.first;
        auto it = sum_by_channel.find(ch);
        if (it == sum_by_channel.end())
            continue;
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
            mc_total_.reset(static_cast<TH1D*>(uptr->Clone((spec_.id + "_mc_total").c_str())));
            mc_total_->SetDirectory(nullptr);
        } else {
            mc_total_->Add(uptr.get());
        }
    }

    if (!data_.empty()) {
        std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
        for (size_t ie = 0; ie < data_.size(); ++ie) {
            const Entry* e = data_[ie];
            if (!e)
                continue;
            auto n0 = selection::apply(e->rnode(), spec_.sel, *e);
            auto n = (spec_.expr.empty() ? n0 : n0.Define("_rx_expr_", spec_.expr));
            const std::string var = spec_.expr.empty() ? spec_.id : "_rx_expr_";
            parts.push_back(n.Histo1D(spec_.model("_data_src" + std::to_string(ie)), var));
        }
        for (auto& rr : parts) {
            const TH1D& h = rr.GetValue();
            if (!data_hist_) {
                data_hist_.reset(static_cast<TH1D*>(h.Clone((spec_.id + "_data").c_str())));
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
        sig->SetLineColor(kGreen + 2);
        sig->SetLineStyle(kDashed);
        sig->SetLineWidth(2);
        sig->SetFillStyle(0);
        sig_hist_ = std::move(sig);
    }
}

void rarexsec::plot::StackedHist::draw_stack_and_unc(TPad* p_main, double& max_y) {
    if (!p_main)
        return;
    p_main->cd();

    if (auto* hists = stack_->GetHists()) {
        for (TObject* obj = hists->First(); obj != nullptr; obj = hists->After(obj)) {
            if (auto* hist = dynamic_cast<TH1*>(obj)) {
                hist->SetLineColor(kBlack);
            }
        }
    }

    stack_->Draw("HIST");
    TH1* frame = stack_->GetHistogram();
    if (frame) {
        frame->SetLineWidth(2);
    }
    if (frame && spec_.xmin < spec_.xmax) {
        frame->GetXaxis()->SetRangeUser(spec_.xmin, spec_.xmax);
    }
    if (frame) {
        frame->GetXaxis()->SetNdivisions(510);
        frame->GetXaxis()->SetTickLength(0.02);
        if (!opt_.x_title.empty())
            frame->GetXaxis()->SetTitle(opt_.x_title.c_str());
        if (!opt_.y_title.empty())
            frame->GetYaxis()->SetTitle(opt_.y_title.c_str());
    }
    if (mc_total_) {
        if (opt_.total_cov || !opt_.syst_bin.empty()) {
            apply_total_errors(*mc_total_, opt_.total_cov.get(),
                               opt_.syst_bin.empty() ? nullptr : &opt_.syst_bin);
        }

        max_y = mc_total_->GetMaximum() + mc_total_->GetBinError(mc_total_->GetMaximumBin());
        if (opt_.y_max > 0)
            max_y = opt_.y_max;

        stack_->SetMaximum(max_y * (opt_.use_log_y ? 10. : 1.3));
        stack_->SetMinimum(opt_.use_log_y ? 0.1 : opt_.y_min);

        auto* h = static_cast<TH1D*>(mc_total_->Clone((spec_.id + "_mc_totband").c_str()));
        h->SetDirectory(nullptr);
        h->SetFillColor(kBlack);
        h->SetFillStyle(3004);
        h->SetMarkerSize(0);
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
        // Draw shaded band plus a crisp outline like the reference
        h->Draw("E2 SAME");
        h->Draw("E1 SAME");
    }
    if (sig_hist_)
        sig_hist_->Draw("HIST SAME");
    if (data_hist_)
        data_hist_->Draw("E1 SAME");
}

void rarexsec::plot::StackedHist::draw_ratio(TPad* p_ratio) {
    if (!p_ratio || !data_hist_ || !mc_total_)
        return;
    p_ratio->cd();

    auto ratio = std::unique_ptr<TH1D>(static_cast<TH1D*>(
        data_hist_->Clone((spec_.id + "_ratio").c_str())));
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
        band.reset(static_cast<TH1D*>(mc_total_->Clone((spec_.id + "_ratio_band").c_str())));
        band->SetDirectory(nullptr);
        const int nb = band->GetNbinsX();
        for (int i = 1; i <= nb; ++i) {
            const double m = mc_total_->GetBinContent(i);
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
    if (!p)
        return;
    p->cd();
    // Use the same box placement as the reference analysis::StackedHistogramPlot
    legend_ = std::make_unique<TLegend>(0.12, 0.0, 0.95, 0.75);
    auto* leg = legend_.get();
    if (!opt_.legend_on_top) {
        leg->SetX1NDC(opt_.leg_x1);
        leg->SetY1NDC(opt_.leg_y1);
        leg->SetX2NDC(opt_.leg_x2);
        leg->SetY2NDC(opt_.leg_y2);
    }
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);

    int n_entries = static_cast<int>(mc_ch_hists_.size());
    if (mc_total_)
        ++n_entries;
    if (sig_hist_)
        ++n_entries;
    if (data_hist_)
        ++n_entries;
    if (n_entries > 0) {
        const int n_cols = (n_entries > 4) ? 3 : 2;
        leg->SetNColumns(n_cols);

        // Scale the text size with the number of rows to maintain a consistent look
        const int n_rows = (n_entries + n_cols - 1) / n_cols;
        const double usable_h = 0.75 * 0.70;
        double txt = usable_h / std::max(1, n_rows);
        txt *= 0.90; // allow for line spacing
        txt = std::clamp(txt, 0.025, 0.050);
        leg->SetTextSize(txt);
    }

    legend_proxies_.clear();

    for (size_t i = 0; i < mc_ch_hists_.size(); ++i) {
        int ch = chan_order_.at(i);
        double sum = mc_ch_hists_[i]->Integral();
        std::string label = rarexsec::plot::Channels::label(ch);
        // Match the reference convention for empty set
        if (label == "#emptyset")
            label = "\xE2\x88\x85";
        if (opt_.annotate_numbers) {
            label += " : " + rarexsec::plot::Plotter::fmt_commas(sum, 2);
        }
        auto proxy = std::unique_ptr<TH1D>(static_cast<TH1D*>(
            mc_ch_hists_[i]->Clone((spec_.id + "_leg_ch" + std::to_string(ch)).c_str())));
        proxy->SetDirectory(nullptr);
        proxy->Reset("ICES");

        auto* entry = leg->AddEntry(proxy.get(), label.c_str(), "f");
        leg->SetEntrySeparation(0.01);
        legend_proxies_.push_back(std::move(proxy));
        (void)entry;
    }

    if (mc_total_) {
        auto proxy = std::unique_ptr<TH1D>(static_cast<TH1D*>(
            mc_total_->Clone((spec_.id + "_leg_unc").c_str())));
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
        // Reference legend uses just "Signal"
        leg->AddEntry(sig_hist_.get(), "Signal", "l");
    }

    if (data_hist_) {
        leg->AddEntry(data_hist_.get(), "Data", "lep");
    }

    leg->Draw();
}

void rarexsec::plot::StackedHist::draw_cuts(TPad* p, double max_y) {
    if (!opt_.show_cuts || opt_.cuts.empty())
        return;
    p->cd();
    TH1* frame = stack_->GetHistogram();
    if (!frame)
        return;
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
    if (!p)
        return;
    p->cd();

    const std::string line1 = "#bf{#muBooNE Simulation, Preliminary}";

    const auto sum_pot = [](const std::vector<const rarexsec::Entry*>& entries) {
        constexpr double rel_tol = 1e-6;
        double total = 0.0;
        std::map<std::pair<std::string, std::string>, std::vector<double>> seen;
        for (const auto* e : entries) {
            if (!e || e->pot_nom <= 0.0)
                continue;
            const auto key = std::make_pair(e->beamline, e->period);
            auto& values = seen[key];
            const double pot = e->pot_nom;
            const bool already_accounted = std::any_of(values.begin(), values.end(), [&](double v) {
                const double scale = std::max(1.0, std::max(std::abs(v), std::abs(pot)));
                return std::abs(v - pot) <= rel_tol * scale;
            });
            if (!already_accounted) {
                values.push_back(pot);
                total += pot;
            }
        }
        return total;
    };

    double pot_value = opt_.total_protons_on_target;
    if (pot_value <= 0.0)
        pot_value = sum_pot(data_);
    if (pot_value <= 0.0)
        pot_value = sum_pot(mc_);

    auto format_pot = [](double value) {
        std::ostringstream ss;
        ss << std::scientific << std::setprecision(2) << value;
        auto text = ss.str();
        const auto pos = text.find('e');
        if (pos != std::string::npos) {
            const int exponent = std::stoi(text.substr(pos + 1));
            text = text.substr(0, pos) + " #times 10^{" + std::to_string(exponent) + "}";
        }
        return text;
    };

    const std::string pot_str = pot_value > 0.0 ? format_pot(pot_value) : "N/A";

    const auto first_non_empty = [](const std::vector<const rarexsec::Entry*>& entries,
                                    auto getter) -> std::string {
        for (const auto* e : entries) {
            if (!e)
                continue;
            auto value = getter(*e);
            if (!value.empty())
                return value;
        }
        return {};
    };

    auto beam_name = opt_.beamline;
    if (beam_name.empty())
        beam_name = first_non_empty(data_, [](const auto& e) { return e.beamline; });
    if (beam_name.empty())
        beam_name = first_non_empty(mc_, [](const auto& e) { return e.beamline; });
    if (beam_name == "numi_fhc")
        beam_name = "NuMI FHC";
    else if (beam_name == "numi_rhc")
        beam_name = "NuMI RHC";
    if (beam_name.empty())
        beam_name = "N/A";

    std::vector<std::string> runs = opt_.run_numbers;
    if (runs.empty())
        runs = opt_.periods;
    if (runs.empty()) {
        auto run = first_non_empty(data_, [](const auto& e) { return e.period; });
        if (run.empty())
            run = first_non_empty(mc_, [](const auto& e) { return e.period; });
        if (!run.empty())
            runs.push_back(std::move(run));
    }

    auto format_run = [](std::string label) {
        if (label.rfind("run", 0) == 0)
            label.erase(0, 3);
        try {
            label = rarexsec::plot::Plotter::fmt_commas(std::stod(label), 0);
        } catch (...) {
        }
        return label;
    };

    std::string runs_str = "N/A";
    if (!runs.empty()) {
        std::ostringstream ss;
        for (size_t i = 0; i < runs.size(); ++i) {
            if (i)
                ss << ", ";
            ss << format_run(runs[i]);
        }
        runs_str = ss.str();
    }

    std::string region_label = opt_.analysis_region_label;
    if (region_label.empty())
        region_label = selection_label(spec_.sel);

    const std::string line2 = "Beam(s), Run(s): " + beam_name + ", " + runs_str +
                              " (" + pot_str + " POT)";
    const std::string line3 = "Analysis Region: " + region_label + " (" +
                              rarexsec::plot::Plotter::fmt_commas(total_mc, 2) + " events)";

    TLatex watermark;
    watermark.SetNDC();
    watermark.SetTextAlign(33);
    watermark.SetTextFont(62);
    watermark.SetTextSize(0.05);
    const double x = 1 - p->GetRightMargin() - 0.03;
    const double top = 1 - p->GetTopMargin();
    watermark.DrawLatex(x, top - 0.03, line1.c_str());
    watermark.SetTextFont(42);
    watermark.SetTextSize(0.05 * 0.8);
    watermark.DrawLatex(x, top - 0.09, line2.c_str());
    watermark.DrawLatex(x, top - 0.15, line3.c_str());
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
    if (want_ratio())
        draw_ratio(p_ratio);
    if (p_main)
        p_main->RedrawAxis();
    canvas.Update();
}

void rarexsec::plot::StackedHist::draw_and_save(const std::string& image_format) {
    std::filesystem::create_directories(output_directory_);
    TCanvas canvas(plot_name_.c_str(), plot_name_.c_str(), 800, 600);
    draw(canvas);
    const std::string fmt = image_format.empty() ? "png" : image_format;
    const std::string out = output_directory_ + "/" + plot_name_ + "." + fmt;
    if (fmt == "pdf") {
        canvas.SaveAs(out.c_str());
        return;
    }

    std::unique_ptr<TImage> image(TImage::Create());
    if (image) {
        image->FromPad(&canvas);
        image->WriteImage(out.c_str());
        return;
    }

    // Fallback to ROOT's SaveAs when TImage is unavailable (e.g. ASImage not
    // built). This ensures we still produce an output file instead of failing
    // silently.
    canvas.SaveAs(out.c_str());
}
