#pragma once

#include <memory>
#include <string>
#include <vector>

#include <TH1.h>

#include "rarexsec/hist/HistogramPlot.hh"

namespace rarexsec::hist {

class StackedHistogramPlot : public HistogramPlot {
public:
    using histogram_ptr = std::unique_ptr<TH1>;

    StackedHistogramPlot(std::string plot_name,
                         std::string output_directory = "plots");

    void add_histogram(histogram_ptr histogram, int color, std::string label);
    void set_axis_labels(std::string x_axis, std::string y_axis);
    void set_legend_title(std::string title);

protected:
    void draw(TCanvas& canvas) override;

private:
    struct Entry {
        histogram_ptr histogram;
        int color;
        std::string label;
    };

    std::vector<Entry> entries_;
    std::string x_axis_title_;
    std::string y_axis_title_;
    std::string legend_title_;
};

} // namespace rarexsec::hist

