#pragma once
#include <memory>
#include <string>
#include <vector>
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"

#include "rarexsec/plot/Plotter.hh"

namespace rarexsec {
namespace plot {

class StackedHist {
public:
    StackedHist(H1Spec spec,
                Options opt,
                std::vector<const Entry*> mc,
                std::vector<const Entry*> data);
    ~StackedHist() = default;

    void draw_and_save(const std::string& image_format);

protected:
    void draw(TCanvas& canvas);

private:
    bool want_ratio() const { return opt_.show_ratio && data_hist_ && mc_total_; }
    void build_histograms();
    void setup_pads(TCanvas& c, TPad*& p_main, TPad*& p_ratio) const;
    void draw_stack_and_unc(TPad* p_main, double& max_y);
    void draw_ratio(TPad* p_ratio);
    void draw_legend(TPad* p);
    void draw_cuts(TPad* p, double max_y);
    void draw_watermark(TPad* p, double total_mc) const;

    H1Spec spec_;
    Options opt_;
    std::vector<const Entry*> mc_;
    std::vector<const Entry*> data_;
    std::string plot_name_;
    std::string output_directory_;
    std::unique_ptr<THStack> stack_;
    std::vector<std::unique_ptr<TH1D>> mc_ch_hists_;
    std::unique_ptr<TH1D> mc_total_;
    std::unique_ptr<TH1D> data_hist_;
    std::unique_ptr<TH1D> sig_hist_;
    std::vector<int> chan_order_;
};

}
}