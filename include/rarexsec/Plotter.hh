#pragma once
#include <cctype>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "rarexsec/proc/Selection.hh"
#include "rarexsec/Hub.hh"

namespace rarexsec {
namespace plot {

enum class CutDir { LessThan, GreaterThan };

struct CutLine {
    double x;
    CutDir dir = CutDir::GreaterThan;
};

struct H1Spec {
    std::string id;
    std::string title;
    std::string expr;
    std::string weight = "w_nominal";
    int nbins = 10;
    double xmin = 0.;
    double xmax = 1.;
    selection::Preset sel = selection::Preset::Empty;

    ROOT::RDF::TH1DModel model(const std::string& suffix) const {
        return ROOT::RDF::TH1DModel((id + suffix).c_str(), title.c_str(), nbins, xmin, xmax);
    }
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
    bool show_cuts = true;
    std::vector<CutLine> cuts;
    std::string beamline;
    std::vector<std::string> periods;
    double leg_x1 = 0.60, leg_y1 = 0.55, leg_x2 = 0.88, leg_y2 = 0.88;
};

class Plotter {
public:
    explicit Plotter(Options opt = {}) : opt_(std::move(opt)) {}
    void draw_stack_by_channel(const H1Spec& spec,
                               const std::vector<const Entry*>& mc,
                               const std::vector<const Entry*>& data = {}) const;

    virtual void set_global_style() const;

    static std::string fmt_commas(double v, int prec = -1) {
        std::ostringstream s;
        if (prec >= 0) s << std::fixed << std::setprecision(prec);
        s << v;
        auto x = s.str();
        auto pos = x.find('.');
        for (int i = (pos == std::string::npos ? static_cast<int>(x.size()) : static_cast<int>(pos)) - 3;
             i > 0; i -= 3) x.insert(i, ",");
        return x;
    }

    static std::string sanitise(std::string s) {
        for (char& c : s) {
            if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '_' || c == '.')) c = '_';
        }
        return s;
    }

private:
    Options opt_;
};

}
}