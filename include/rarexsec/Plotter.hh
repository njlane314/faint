#pragma once
#include <string>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "rarexsec/Selection.hh"
#include "rarexsec/Hub.hh"

namespace rarexsec::plot {

struct Hist1D {
    std::string name;
    std::string title;
    std::string expr;
    std::string weight = "w_nominal";
    int nbins = 10;
    double xmin = 0.;
    double xmax = 1.;
    selection::Preset sel = selection::Preset::Empty;
};

enum class CutDir { LessThan, GreaterThan };

struct CutLine {
    double x;
    CutDir dir = CutDir::GreaterThan;
};

struct Options {
    std::string out_dir = "plots";
    std::string image_format = "png";
    bool show_ratio = true;
    bool use_log_y = false;
    bool annotate_numbers = true;
    bool overlay_signal = false;
    std::vector<int> signal_channels = {15,16};
    double y_min = 0.;
    double y_max = -1.;
    std::vector<CutLine> cuts;
    std::string beamline;
    std::vector<std::string> periods;
};

class Plotter {
public:
    explicit Plotter(Options opt = {}) : opt_(std::move(opt)) {}
    void draw_stack_by_channel(const Hist1D& spec,
                               const std::vector<const Entry*>& mc,
                               const std::vector<const Entry*>& data = {}) const;
private:
    Options opt_;
};

}
