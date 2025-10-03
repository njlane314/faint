#include "rarexsec/plot/Plotter.hh"
#include "rarexsec/plot/StackedByChannelPlot.hh"
#include "rarexsec/Channels.hh"

namespace rarexsec::plot {

void Plotter::draw_stack_by_channel(const Hist1D& spec,
                                    const std::vector<const Entry*>& mc,
                                    const std::vector<const Entry*>& data) const {
    auto order = rarexsec::Channels::order();
    StackedByChannelPlot plot(spec.name, opt_.out_dir, spec, opt_, mc, data, order);
    plot.drawAndSave(opt_.image_format);
}

}
