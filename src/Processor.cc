#include "rarexsec/Processor.hh"
#include "rarexsec/Hub.hh"
#include <cmath>

ROOT::RDF::RNode rarexsec::Processor::run(ROOT::RDF::RNode node, const rarexsec::Entry& rec) const {
    const bool is_data = (rec.kind == rarexsec::sample::origin::data);
    const bool is_ext = (rec.kind == rarexsec::sample::origin::ext);
    const bool is_mc = !is_data && !is_ext;

    node = node.Define("is_data",         [is_data]{ return is_data; });
    node = node.Define("is_mc",           [is_mc]{ return is_mc; });
    node = node.Define("is_ext",          [is_ext]{ return is_ext; }); 

    double scale_mc  = 1.0;
    if (is_mc && rec.pot > 0.0 && rec.pot_eff > 0.0) 
        scale_mc = rec.pot_eff / rec.pot;                  

    double scale_ext = 1.0;
    if (is_ext && rec.trig > 0.0 && rec.trig_eff > 0.0) 
        scale_ext = rec.trig_eff / rec.trig;            

    node = node.Define("w_base", [is_mc, is_ext, scale_mc, scale_ext] {
        return is_mc ? scale_mc : (is_ext ? scale_ext : 1.0); 
    });

    if (is_mc) {
        node = node.Define(
            "w_nominal",
            [](double w, double w_spline, double w_tune) {
                double out = w * w_spline * w_tune;
                if (!std::isfinite(out) || out < 0.0) return 1.0;
                return out;
            },
            {"w_base", "weightSpline", "weightTune"});
    } else {
        node = node.Define("w_nominal", [](double w) { return w; }, {"w_base"});
    }

    if (is_mc) {
        node = node.Define(
            "in_fiducial",
            [](float x, float y, float z) {
                return rarexsec::fiducial::is_in_truth_volume(x, y, z);
            },
            {"neutrino_vertex_x", "neutrino_vertex_y", "neutrino_vertex_z"});

        node = node.Define(
            "is_strange",
            [](int kplus, int kminus, int kzero, int lambda0, int sigplus, int sigzero, int sigminus) {
                return (kplus + kminus + kzero + lambda0 + sigplus + sigzero + sigminus) > 0;
            },
            {"count_kaon_plus", "count_kaon_minus", "count_kaon_zero", "count_lambda", "count_sigma_plus", "count_sigma_zero", "count_sigma_minus"});
    } 

    return node;
}

const rarexsec::Processor& rarexsec::processor() {
    static const Processor ep{};
    return ep;
}
