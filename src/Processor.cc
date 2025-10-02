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
            
        node = node.Define(
            "interation_mode",
            [](int mode) {
                switch (mode) {
                    case 0:  return 0;   // QEL
                    case 1:  return 1;   // RES
                    case 2:  return 2;   // DIS
                    case 3:  return 3;   // COH
                    case 10: return 10;  // MEC/2p2h
                    default: return -1;  // uncategorised
                }
            },
            {"interaction_mode"});

        node = node.Define(
            "analysis_channels",
            [](bool fv, int nu, int ccnc, int s, int npi, int np, int npi0, int ngamma) {
                if (!fv) {
                    if (nu == 0) return 1;   // out-of-fv, no neutrino?
                    return 2;                // out-of-fv, with neutrino
                }
                if (ccnc == 1) return 14;    // NC in-fv
                if (ccnc == 0 && s > 0) {
                    if (s == 1) return 15;   // CC + 1 strange
                    return 16;               // CC + >1 strange
                }
                if (std::abs(nu) == 12 && ccnc == 0) return 17; // ν_e CC (no strange)
                if (std::abs(nu) == 14 && ccnc == 0) {          // ν_μ CC (no strange)
                    if (npi == 0 && np > 0)       return 10;    // CC0π, ≥1p
                    if (npi == 1 && npi0 == 0)    return 11;    // CC1π± only
                    if (npi0 > 0 || ngamma >= 2)  return 12;    // CC π0 / γ-rich
                    if (npi > 1)                  return 13;    // CC Nπ± (N>1)
                    return 18;                                  // CC other
                }
                return 99; // other
            },
            {"in_fiducial", "neutrino_pdg", "interaction_ccnc", "mc_n_strange", "mc_n_pion", "mc_n_proton", "count_pi_zero", "count_gamma"});

            node = node.Define(
                "is_signal",
                [](int ch){ return ch == 15 || ch == 16; },
                {"channel_definitions"});

            node = node.Define(
                "recognised_signal",
                [](bool is_sig, float purity, float completeness) {
                    return is_sig && purity > 0.5f && completeness > 0.1f;
                },
                {"is_signal", "neutrino_purity_from_pfp", "neutrino_completeness_from_pfp"});
    } 

    return node;
}

const rarexsec::Processor& rarexsec::processor() {
    static const Processor ep{};
    return ep;
}
