#include "rarexsec/hist/StackedHistogramPlot.hh"

#include <cstddef>
#include <stdexcept>

#include <THStack.h>
#include <TLegend.h>

namespace rarexsec::hist {

StackedHistogramPlot::StackedHistogramPlot(std::string plot_name,
                                           std::string output_directory)
    : HistogramPlot(std::move(plot_name), std::move(output_directory)) {}

void StackedHistogramPlot::add_histogram(histogram_ptr histogram, int color, std::string label) {
    if (!histogram) {
        throw std::invalid_argument("Histogram pointer cannot be null");
    }

    histogram->SetDirectory(nullptr);
    histogram->SetLineColor(color);
    histogram->SetMarkerColor(color);
    histogram->SetFillColor(color);
    histogram->SetFillStyle(1001);

    entries_.push_back(Entry{std::move(histogram), color, std::move(label)});
}

void StackedHistogramPlot::set_axis_labels(std::string x_axis, std::string y_axis) {
    x_axis_title_ = std::move(x_axis);
    y_axis_title_ = std::move(y_axis);
}

void StackedHistogramPlot::set_legend_title(std::string title) {
    legend_title_ = std::move(title);
}

void StackedHistogramPlot::draw(TCanvas& canvas) {
    if (entries_.empty()) {
        throw std::runtime_error("No histograms provided for stacked plot '" + name() + "'");
    }

    canvas.cd();
    THStack stack(("stack_" + sanitise(name())).c_str(), name().c_str());

    TLegend legend(0.62, 0.60, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    if (!legend_title_.empty()) {
        legend.SetHeader(legend_title_.c_str(), "C");
    }

    std::size_t index = 0;
    for (auto& entry : entries_) {
        auto* histogram = entry.histogram.get();
        if (!histogram) {
            continue;
        }
        histogram->SetLineColor(entry.color);
        histogram->SetMarkerColor(entry.color);
        histogram->SetFillColor(entry.color);
        histogram->SetFillStyle(1001);
        histogram->SetName(("component_" + std::to_string(index++) + "_" + sanitise(entry.label)).c_str());
        stack.Add(histogram);
        legend.AddEntry(histogram, entry.label.c_str(), "f");
    }

    stack.Draw("hist");
    if (!x_axis_title_.empty()) {
        stack.GetXaxis()->SetTitle(x_axis_title_.c_str());
    }
    if (!y_axis_title_.empty()) {
        stack.GetYaxis()->SetTitle(y_axis_title_.c_str());
    }

    legend.Draw();
}

} // namespace rarexsec::hist

