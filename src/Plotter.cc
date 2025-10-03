#include "rarexsec/Plotter.hh"
#include "rarexsec/plot/StackedHist.hh"

#include "TROOT.h"
#include "TStyle.h"

void rarexsec::plot::Plotter::draw_stack_by_channel(const H1Spec& spec,
                                    const std::vector<const Entry*>& mc,
                                    const std::vector<const Entry*>& data) const {
    set_global_style();
    StackedHist plot(spec, opt_, mc, data);
    plot.draw_and_save(opt_.image_format);
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
    style->SetPadTickX(1);
    style->SetPadTickY(1);
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
    style->SetTitleFillColor(0);
    style->SetTitleBorderSize(0);
    gROOT->SetStyle("PlotterStyle");
    gROOT->ForceStyle();
}