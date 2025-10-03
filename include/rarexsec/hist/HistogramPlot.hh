#pragma once

#include <memory>
#include <string>

#include <TCanvas.h>
#include <TImage.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

namespace rarexsec::hist {

class HistogramPlot {
public:
    explicit HistogramPlot(std::string plot_name,
                           std::string output_directory = "plots")
        : plot_name_(std::move(plot_name)),
          output_directory_(std::move(output_directory)) {
        gSystem->mkdir(output_directory_.c_str(), true);
    }

    virtual ~HistogramPlot() = default;

    void draw_and_save(const std::string& format = "png") {
        gROOT->SetBatch(kTRUE);
        gSystem->mkdir(output_directory_.c_str(), true);
        set_global_style();

        TCanvas canvas(plot_name_.c_str(), plot_name_.c_str(), 800, 600);
        draw(canvas);
        canvas.Update();

        const auto output_path = output_directory_ + "/" + plot_name_ + "." + format;
        if (format == "pdf") {
            canvas.SaveAs(output_path.c_str());
            return;
        }

        std::unique_ptr<TImage> image{TImage::Create()};
        if (image) {
            image->FromPad(&canvas);
            image->WriteImage(output_path.c_str());
        }
    }

    static std::string sanitise(std::string value) {
        for (auto& c : value) {
            if (c == '.' || c == '/' || c == ' ') {
                c = '_';
            }
        }
        return value;
    }

protected:
    [[nodiscard]] const std::string& name() const noexcept { return plot_name_; }
    [[nodiscard]] const std::string& output_directory() const noexcept { return output_directory_; }

    virtual void draw(TCanvas& canvas) = 0;

    virtual void set_global_style() const {
        constexpr int font_style = 42;
        if (!gROOT->GetStyle("RarexsecPlot")) {
            auto style = std::make_unique<TStyle>("RarexsecPlot", "Rarexsec Plot Style");
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
            gROOT->GetListOfStyles()->Add(style.release());
        }

        gROOT->SetStyle("RarexsecPlot");
        gROOT->ForceStyle();
    }

private:
    std::string plot_name_;
    std::string output_directory_;
};

} // namespace rarexsec::hist

