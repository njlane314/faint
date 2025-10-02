#include "Processor.hh"

ROOT::RDF::RNode rarexsec::Processor::run(ROOT::RDF::RNode node,
                                          rarexsec::sample::origin kind) const {
    const bool is_data = (kind == rarexsec::sample::origin::data);

    node = node.Define("is_data",         [is_data]{ return is_data; });
    node = node.Define("is_simulation",   [is_data]{ return !is_data; });

    node = node.Define("w_nominal", [is_data]{ return is_data ? 1.0 : 1.0; });

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