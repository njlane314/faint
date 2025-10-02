#include "Processor.hh"

ROOT::RDF::RNode rarexsec::Processor::run(ROOT::RDF::RNode node, const rarexsec::Entry rec) const {
    const bool is_data = (rec.kind == rarexsec::sample::origin::data);
    const bool is_ext = (rec.kind == rarexsec::sample::origin::ext);

    node = node.Define("is_data",         [is_data]{ return is_data; });
    node = node.Define("is_simulation",   [is_data]{ return !is_data; });

    node = node.Define("w_nominal", []{ return 1.0; }); 

    if (!is_data || !is_ext) {
        double scale = 1.0;
        if (rec.pot > 0.0 && rec.pot_eff > 0.0)
            scale = rec.pot_eff / rec.pot;
        node = node.Define("w_base", [scale]() { return is_data ? scale : 1.0; });
    
        node = node.Define("w_nominal",
            [](double w, float w_spline, float w_tune) {
                double out = w;
                out *= w_spline * w_tune;
                
                if (!std::isfinite(out) || out < 0)
                    return 1.0;
                return out;
            },
        {"w_base", "weightSpline", "weightTune"});
    } else if (is_ext)
        double scale = 1.0;
        if (rec.trigs > 0 && rec.trigs_eff > 0) 
            scale = static_cast<double>(rec.trigs_eff) / static_cast<double>(sample_triggers_);
        node = node.Define("w_nominal", [scale]() { return scale; });
    }


    // -------------------------
    // YOUR ANALYSIS LOGIC HERE
    // e.g.,
    // if (!is_data) {
    //   node = node.Define("w_genie", "genie_w");      // if such a column exists
    //   node = node.Define("w", "w_nominal * w_genie");
    // } else {
    //   node = node.Define("w", "w_nominal");
    // }
    // node = node.Filter("n_tracks > 0", ">=1 reco track");
    // -------------------------

    return node;
}

const Processor& processor() {
    static const Processor ep{};
    return ep;
}

} 