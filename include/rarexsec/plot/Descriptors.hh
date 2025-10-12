#pragma once

#include <cctype>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RDFHelpers.hxx>
#include <TMatrixDSym.h>

#include "rarexsec/proc/Selection.hh"

namespace rarexsec {
namespace plot {

enum class CutDir { LessThan,
                    GreaterThan };

struct CutSpec {
    double x = 0.0;
    CutDir dir = CutDir::GreaterThan;
};

struct Options {
    std::string out_dir = ".";
    std::string image_format = "png";
    bool show_ratio = false;
    bool show_ratio_band = false;
    bool use_log_y = false;
    bool annotate_numbers = false;
    bool overlay_signal = false;
    bool legend_on_top = false;
    bool show_cuts = false;
    double legend_split = 0.75;
    double y_min = 0.0;
    double y_max = -1.0;
    std::string x_title;
    std::string y_title;
    std::vector<int> signal_channels;
    std::shared_ptr<TMatrixDSym> total_cov;
    std::vector<double> syst_bin;
    std::vector<CutSpec> cuts;
    double total_protons_on_target = 0.0;
    std::string beamline;
    std::vector<std::string> run_numbers;
    std::vector<std::string> periods;
    std::string analysis_region_label;
};

struct TH1DModel {
    std::string id;
    std::string name;
    std::string title;
    std::string expr;
    std::string weight = "w_nominal";
    int nbins = 1;
    double xmin = 0.0;
    double xmax = 1.0;
    selection::Preset sel = selection::Preset::InclusiveMuCC;

    ROOT::RDF::TH1DModel model(const std::string& suffix = "") const {
        const std::string base = !id.empty() ? id : name;
        const std::string hist_name = sanitise(base + suffix);
        const std::string hist_title = title.empty() ? base : title;
        return ROOT::RDF::TH1DModel(hist_name.c_str(), hist_title.c_str(), nbins, xmin, xmax);
    }

    std::string axis_title() const {
        if (!title.empty()) {
            return title;
        }
        const std::string base = !name.empty() ? name : id;
        if (base.empty()) {
            return ";x;Events";
        }
        return ";" + base + ";Events";
    }

  private:
    static std::string sanitise(const std::string& raw) {
        std::string out;
        out.reserve(raw.size());
        for (unsigned char c : raw) {
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
};

using Histogram1DSpec = TH1DModel;

}
}
