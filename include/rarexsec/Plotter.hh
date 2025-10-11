#pragma once

#include <string>
#include <vector>

#include "TMatrixDSym.h"
#include "rarexsec/Hub.hh"
#include "rarexsec/plot/Descriptors.hh"
#include "rarexsec/plot/EventDisplay.hh"

namespace rarexsec::plot {

class Plotter {
  public:
    Plotter();
    explicit Plotter(Options opt);

    const Options& options() const noexcept;
    Options& options() noexcept;

    void set_options(Options opt);

    void draw_stack_by_channel(const Histogram1DSpec& spec,
                               const std::vector<const Entry*>& mc) const;

    void draw_stack_by_channel(const Histogram1DSpec& spec,
                               const std::vector<const Entry*>& mc,
                               const std::vector<const Entry*>& data) const;

    void draw_unstacked_by_channel(const Histogram1DSpec& spec,
                                   const std::vector<const Entry*>& mc,
                                   bool normalize_to_pdf = true,
                                   int line_width = 3) const;

    void draw_unstacked_by_channel(const Histogram1DSpec& spec,
                                   const std::vector<const Entry*>& mc,
                                   const std::vector<const Entry*>& data,
                                   bool normalize_to_pdf,
                                   int line_width) const;

    void draw_stack_by_channel_with_cov(const Histogram1DSpec& spec,
                                        const std::vector<const Entry*>& mc,
                                        const std::vector<const Entry*>& data,
                                        const TMatrixDSym& total_cov) const;

    void draw_event_display(EventDisplay::Spec spec,
                            EventDisplay::Options opt,
                            EventDisplay::DetectorData data) const;

    void draw_event_display(EventDisplay::Spec spec,
                            EventDisplay::Options opt,
                            EventDisplay::SemanticData data) const;

    static std::string sanitise(const std::string& name);
    static std::string fmt_commas(double value, int precision);

    void set_global_style() const;

  private:
    Options opt_;
};

}
