#include "rarexsec/Plotter.hh"
#include "rarexsec/plot/StackedHist.hh"

void rarexsec::plot::Plotter::draw_stack_by_channel(const H1Spec& spec,
                                    const std::vector<const Entry*>& mc,
                                    const std::vector<const Entry*>& data) const {
    StackedHist plot(spec, opt_, mc, data);
    plot.draw_and_save(opt_.image_format);
}