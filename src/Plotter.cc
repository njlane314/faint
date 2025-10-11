#include "rarexsec/Plotter.hh"
#include <cctype>
#include <iomanip>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>
#include <TROOT.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TMatrixDSym.h>
#include "rarexsec/plot/EventDisplay.hh"
#include "rarexsec/plot/StackedHist.hh"
#include "rarexsec/plot/UnstackedHist.hh"

rarexsec::plot::Plotter::Plotter() = default;

rarexsec::plot::Plotter::Plotter(rarexsec::plot::Options opt)
    : opt_(std::move(opt)) {}

const rarexsec::plot::Options& rarexsec::plot::Plotter::options() const noexcept { return opt_; }

rarexsec::plot::Options& rarexsec::plot::Plotter::options() noexcept { return opt_; }

void rarexsec::plot::Plotter::set_options(rarexsec::plot::Options opt) { opt_ = std::move(opt); }

void Plotter::set_options(Options opt) { opt_ = std::move(opt); }

void Plotter::draw_stack_by_channel(const Histogram1DSpec& spec,
                                    const std::vector<const Entry*>& mc) const {
    static const std::vector<const Entry*> empty_data{};
    draw_stack_by_channel(spec, mc, empty_data);
}

void Plotter::draw_stack_by_channel(const Histogram1DSpec& spec,
                                    const std::vector<const Entry*>& mc,
                                    const std::vector<const Entry*>& data) const {
    set_global_style();
    rarexsec::plot::StackedHist plot(spec, opt_, mc, data);
    plot.draw_and_save(opt_.image_format);
}

void Plotter::draw_unstacked_by_channel(const Histogram1DSpec& spec,
                                        const std::vector<const Entry*>& mc,
                                        bool normalize_to_pdf,
                                        int line_width) const {
    static const std::vector<const Entry*> empty_data{};
    draw_unstacked_by_channel(spec, mc, empty_data, normalize_to_pdf, line_width);
}

void Plotter::draw_unstacked_by_channel(const Histogram1DSpec& spec,
                                        const std::vector<const Entry*>& mc,
                                        const std::vector<const Entry*>& data,
                                        bool normalize_to_pdf,
                                        int line_width) const {
    set_global_style();
    rarexsec::plot::UnstackedHist plot(spec, opt_, mc, data, normalize_to_pdf, line_width);
    plot.draw_and_save(opt_.image_format);
}

void Plotter::draw_stack_by_channel_with_cov(const Histogram1DSpec& spec,
                                             const std::vector<const Entry*>& mc,
                                             const std::vector<const Entry*>& data,
                                             const TMatrixDSym& total_cov) const {
    set_global_style();
    auto opt2 = opt_;
    opt2.total_cov = std::make_shared<TMatrixDSym>(total_cov);
    rarexsec::plot::StackedHist plot(spec, std::move(opt2), mc, data);
    plot.draw_and_save(opt2.image_format);
}

void rarexsec::plot::Plotter::draw_event_display(rarexsec::plot::EventDisplay::Spec spec,
                                                 rarexsec::plot::EventDisplay::Options opt,
                                                 rarexsec::plot::EventDisplay::DetectorData data) const {
    rarexsec::plot::EventDisplay display(std::move(spec), std::move(opt), std::move(data));
    display.draw_and_save(opt_.image_format);
}

void rarexsec::plot::Plotter::draw_event_display(rarexsec::plot::EventDisplay::Spec spec,
                                                 rarexsec::plot::EventDisplay::Options opt,
                                                 rarexsec::plot::EventDisplay::SemanticData data) const {
    rarexsec::plot::EventDisplay display(std::move(spec), std::move(opt), std::move(data));
    display.draw_and_save(opt_.image_format);
}

std::string Plotter::sanitise(const std::string& name) {
    std::string out;
    out.reserve(name.size());
    for (unsigned char c : name) {
        if (std::isalnum(c) || c == '_' || c == '-') {
            out.push_back(static_cast<char>(c));
        } else if (c == '.' || c == ' ') {
            out.push_back('_');
        } else {
            out.push_back('_');
        }
    }
    if (out.empty()) {
        return "plot";
    }
    return out;
}

std::string rarexsec::plot::Plotter::fmt_commas(double value, int precision) {
    std::ostringstream ss;
    if (precision >= 0) {
        ss << std::fixed << std::setprecision(precision);
    }
    ss << value;
    std::string text = ss.str();
    const auto pos = text.find('.');
    std::string integer = pos == std::string::npos ? text : text.substr(0, pos);
    std::string fraction = pos == std::string::npos ? std::string{} : text.substr(pos);
    bool negative = false;
    if (!integer.empty() && integer.front() == '-') {
        negative = true;
        integer.erase(integer.begin());
    }
    std::string with_commas;
    for (std::size_t i = 0; i < integer.size(); ++i) {
        if (i != 0 && (integer.size() - i) % 3 == 0) {
            with_commas.push_back(',');
        }
        with_commas.push_back(integer[i]);
    }
    if (negative) {
        with_commas.insert(with_commas.begin(), '-');
    }
    return with_commas + fraction;
}

void rarexsec::plot::Plotter::set_global_style() const {
    const int font_style = 42;
    TStyle* style = new TStyle("PlotterStyle", "Plotter Style");
    style->SetTitleFont(font_style, "X");
    style->SetTitleFont(font_style, "Y");
    style->SetTitleFont(font_style, "Z");
    style->SetTitleSize(0.05, "X");
    style->SetTitleSize(0.05, "Y");
    style->SetTitleSize(0.04, "Z");
    style->SetLabelFont(font_style, "X");
    style->SetLabelFont(font_style, "Y");
    style->SetLabelFont(font_style, "Z");
    style->SetLabelSize(0.045, "X");
    style->SetLabelSize(0.045, "Y");
    style->SetLabelSize(0.045, "Z");
    style->SetTitleOffset(0.93, "X");
    style->SetTitleOffset(1.06, "Y");
    style->SetOptStat(0);
    style->SetOptTitle(0);
    style->SetPadTickX(1);
    style->SetPadTickY(1);
    TGaxis::SetMaxDigits(4);
    style->SetPadLeftMargin(0.15);
    style->SetPadRightMargin(0.05);
    style->SetPadTopMargin(0.07);
    style->SetPadBottomMargin(0.12);
    style->SetMarkerSize(1.0);
    style->SetCanvasColor(0);
    style->SetPadColor(0);
    style->SetFrameFillColor(0);
    style->SetCanvasBorderMode(0);
    style->SetPadBorderMode(0);
    style->SetStatColor(0);
    style->SetFrameBorderMode(0);
    style->SetLineWidth(2);
    style->SetFrameLineWidth(2);
    style->SetTitleFillColor(0);
    style->SetTitleBorderSize(0);
    gROOT->SetStyle("PlotterStyle");
    gROOT->ForceStyle();
}
