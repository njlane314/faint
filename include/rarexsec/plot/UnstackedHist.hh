#pragma once
#include <memory>
#include <string>
#include <vector>

#include "rarexsec/plot/Descriptors.hh"

class TCanvas;
class TLegend;
class TPad;
class TH1D;

namespace rarexsec {
struct Entry;
namespace plot {

class UnstackedHist {
public:
    UnstackedHist(Histogram1DSpec spec,
                  Options opt,
                  std::vector<const Entry*> mc,
                  std::vector<const Entry*> data = {},
                  bool normalize_to_pdf = true,
                  int line_width = 3);

    void draw(TCanvas& canvas);
    void draw_and_save(const std::string& image_format = "");

private:
    void build_histograms();
    void setup_pads(TCanvas& c, TPad*& p_main, TPad*& p_legend) const;
    void draw_curves(TPad* p_main, double& max_y);
    void draw_legend(TPad* p_legend);
    void draw_watermark(TPad* p_main) const;

    Histogram1DSpec spec_;
    Options opt_;
    std::vector<const Entry*> mc_;
    std::vector<const Entry*> data_;
    bool normalize_to_pdf_;
    int line_width_;

    std::string plot_name_;
    std::string output_directory_;

    std::vector<int> chan_order_;
    std::vector<std::unique_ptr<TH1D>> mc_ch_hists_;
    std::unique_ptr<TH1D> data_hist_;

    std::unique_ptr<TLegend> legend_;
    std::vector<std::unique_ptr<TH1D>> legend_proxies_;
};

}
}
