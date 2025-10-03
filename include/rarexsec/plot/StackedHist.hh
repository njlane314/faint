#pragma once
#include <memory>
#include <string>
#include <vector>
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPad.h"
#include "TFile.h"
#include "TLatex.h"
#include "rarexsec/plot/Plotter.hh"
#include "rarexsec/Channels.hh"
#include "rarexsec/Selection.hh"
#include "rarexsec/Hub.hh"

class TCanvas;

namespace rarexsec {
namespace plot {

class StackedHist {
public:
    StackedHist(std::string plot_name,
                         std::string out_dir,
                         Hist1D spec,
                         Options opt,
                         std::vector<const Entry*> mc,
                         std::vector<const Entry*> data,
                         std::vector<int> channel_order);
    ~StackedHist() = default;

    void draw_and_save(const std::string& image_format);

protected:
    void draw(TCanvas& canvas);

private:
    void build_histograms();
    void setup_pads_ratio(TCanvas& c, TPad*& p_main, TPad*& p_ratio) const;
    void setup_pads_legend_top(TCanvas& c, TPad*& p_main, TPad*& p_leg) const;
    void draw_stack_and_unc(TPad* p_main, double& max_y);
    void draw_ratio(TPad* p_ratio);
    void draw_legend(TPad* p, bool in_main);
    void draw_cuts(TPad* p, double max_y);
    void draw_watermark(TPad* p, double total_mc) const;

    Hist1D spec_;
    Options opt_;
    std::vector<const Entry*> mc_;
    std::vector<const Entry*> data_;
    std::vector<int> chan_order_;
    std::string plot_name_;
    std::string output_directory_;
    std::unique_ptr<THStack> stack_;
    std::vector<std::unique_ptr<TH1D>> mc_ch_hists_;
    std::unique_ptr<TH1D> mc_total_;
    std::unique_ptr<TH1D> data_hist_;
    std::unique_ptr<TH1D> sig_hist_;
};

}
}