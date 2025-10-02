#include "rarexsec/Processor.hh"
#include "rarexsec/Hub.hh"
#include "rarexsec/FiducialVolume.hh"
#include "rarexsec/Selection.h"
#include <ROOT/RVec.hxx>
#include <cmath>

ROOT::RDF::RNode rarexsec::Processor::run(ROOT::RDF::RNode node, const rarexsec::Entry& rec) const {
    const bool is_data = (rec.kind == rarexsec::sample::origin::data);
    const bool is_ext = (rec.kind == rarexsec::sample::origin::ext);
    const bool is_mc = !is_data && !is_ext;

    node = node.Define("is_data",         [is_data]{ return is_data; });
    node = node.Define("is_mc",           [is_mc]{ return is_mc; });
    node = node.Define("is_ext",          [is_ext]{ return is_ext; }); 

    double scale_mc  = 1.0;
    if (is_mc && rec.pot_nom > 0.0 && rec.pot_eqv > 0.0)
        scale_mc = rec.pot_eqv / rec.pot_nom;

    double scale_ext = 1.0;
    if (is_ext && rec.trig_nom > 0.0 && rec.trig_eqv > 0.0)
        scale_ext = rec.trig_eqv / rec.trig_nom;

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
            "scattering_mode",
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
            {"analysis_channels"});

        node = node.Define(
            "recognised_signal",
            [](bool is_sig, float purity, float completeness) {
                return is_sig && purity > 0.5f && completeness > 0.1f;
            },
            {"is_signal", "neutrino_purity_from_pfp", "neutrino_completeness_from_pfp"});

    } else {
        const int nonmc_channel = is_data ? 0 : (is_ext ? 1 : 99);
        node = node.Define("in_fiducial", [] { return false; });
        node = node.Define("is_strange", [] { return false; });
        node = node.Define("scattering_mode", [] { return -1; });
        node = node.Define("analysis_channels", [nonmc_channel] { return nonmc_channel; });
        node = node.Define("is_signal", [] { return false; });
        node = node.Define("recognised_signal", [] { return false; });
    }

    node = node.Define(
        "in_reco_fiducial",
        [](float x, float y, float z) {
            return rarexsec::fiducial::is_in_reco_volume(x, y, z);
        },
        {"reco_neutrino_vertex_sce_x", "reco_neutrino_vertex_sce_y", "reco_neutrino_vertex_sce_z"});

    if (!node.HasColumn("n_pfps_gen2")) {
        node = node.Define(
            "n_pfps_gen2",
            [](const ROOT::RVec<unsigned>& gens) {
                return ROOT::VecOps::Sum(gens == 2u);
            },
            {"pfp_generations"});
    }

    if (!node.HasColumn("n_pfps_gen3")) {
        node = node.Define(
            "n_pfps_gen3",
            [](const ROOT::RVec<unsigned>& gens) {
                return ROOT::VecOps::Sum(gens == 3u);
            },
            {"pfp_generations"});
    }

    if (node.HasColumn("track_shower_scores")) {
        node = node.Define(
            "muon_mask",
            [](const ROOT::RVec<float>& scores, const ROOT::RVec<float>& llr,
               const ROOT::RVec<float>& lengths, const ROOT::RVec<float>& dists,
               const ROOT::RVec<float>& start_x, const ROOT::RVec<float>& start_y,
               const ROOT::RVec<float>& start_z, const ROOT::RVec<float>& end_x,
               const ROOT::RVec<float>& end_y, const ROOT::RVec<float>& end_z,
               const ROOT::RVec<unsigned>& gens) {
                ROOT::RVec<bool> mask(scores.size());
                for (std::size_t i = 0; i < scores.size(); ++i) {
                    const bool fid_start = rarexsec::fiducial::is_in_reco_volume(
                        start_x[i], start_y[i], start_z[i]);
                    const bool fid_end = rarexsec::fiducial::is_in_reco_volume(
                        end_x[i], end_y[i], end_z[i]);
                    mask[i] = rarexsec::selection::passes_muon_track_selection(
                        scores[i], llr[i], lengths[i], dists[i], gens[i], fid_start, fid_end);
                }
                return mask;
            },
            {"track_shower_scores", "trk_llr_pid_v", "track_length",
             "track_distance_to_vertex", "track_start_x", "track_start_y",
             "track_start_z", "track_end_x", "track_end_y", "track_end_z",
             "pfp_generations"});

        const auto filter_float = [](const ROOT::RVec<float>& vals,
                                     const ROOT::RVec<bool>& mask) {
            ROOT::RVec<float> out;
            out.reserve(vals.size());
            for (std::size_t i = 0; i < vals.size(); ++i) {
                if (mask[i]) out.push_back(vals[i]);
            }
            return out;
        };

        const auto filter_uint = [](const ROOT::RVec<unsigned>& vals,
                                    const ROOT::RVec<bool>& mask) {
            ROOT::RVec<unsigned> out;
            out.reserve(vals.size());
            for (std::size_t i = 0; i < vals.size(); ++i) {
                if (mask[i]) out.push_back(vals[i]);
            }
            return out;
        };

        const auto filter_costheta = [](const ROOT::RVec<float>& theta,
                                        const ROOT::RVec<bool>& mask) {
            ROOT::RVec<float> out;
            out.reserve(theta.size());
            for (std::size_t i = 0; i < theta.size(); ++i) {
                if (mask[i]) out.push_back(std::cos(theta[i]));
            }
            return out;
        };

        node = node.Define("muon_trk_score_v", filter_float,
                           {"track_shower_scores", "muon_mask"});
        node = node.Define("muon_trk_llr_pid_v", filter_float,
                           {"trk_llr_pid_v", "muon_mask"});
        node = node.Define("muon_trk_start_x_v", filter_float,
                           {"track_start_x", "muon_mask"});
        node = node.Define("muon_trk_start_y_v", filter_float,
                           {"track_start_y", "muon_mask"});
        node = node.Define("muon_trk_start_z_v", filter_float,
                           {"track_start_z", "muon_mask"});
        node = node.Define("muon_trk_end_x_v", filter_float,
                           {"track_end_x", "muon_mask"});
        node = node.Define("muon_trk_end_y_v", filter_float,
                           {"track_end_y", "muon_mask"});
        node = node.Define("muon_trk_end_z_v", filter_float,
                           {"track_end_z", "muon_mask"});
        node = node.Define("muon_trk_length_v", filter_float,
                           {"track_length", "muon_mask"});
        node = node.Define("muon_trk_distance_v", filter_float,
                           {"track_distance_to_vertex", "muon_mask"});
        node = node.Define("muon_pfp_generation_v", filter_uint,
                           {"pfp_generations", "muon_mask"});
        node = node.Define("muon_track_costheta", filter_costheta,
                           {"track_theta", "muon_mask"});

        if (node.HasColumn("n_muons_tot")) {
            node = node.Redefine(
                "n_muons_tot",
                [](const ROOT::RVec<bool>& mask) {
                    return static_cast<unsigned long>(ROOT::VecOps::Sum(mask));
                },
                {"muon_mask"});
        } else {
            node = node.Define(
                "n_muons_tot",
                [](const ROOT::RVec<bool>& mask) {
                    return static_cast<unsigned long>(ROOT::VecOps::Sum(mask));
                },
                {"muon_mask"});
        }

        if (node.HasColumn("has_muon")) {
            node = node.Redefine(
                "has_muon",
                [](unsigned long n_muons) { return n_muons > 0UL; },
                {"n_muons_tot"});
        } else {
            node = node.Define(
                "has_muon",
                [](unsigned long n_muons) { return n_muons > 0UL; },
                {"n_muons_tot"});
        }
    } else {
        if (!node.HasColumn("n_muons_tot")) {
            node = node.Define("n_muons_tot", [] { return 0UL; });
        }
        if (!node.HasColumn("has_muon")) {
            node = node.Define("has_muon", [] { return false; });
        }
    }

    if (node.HasColumn("software_trigger_pre_ext")) {
        node = node.Define(
            "software_trigger",
            [](unsigned run, int pre, int post) {
                return run < 16880 ? pre > 0 : post > 0;
            },
            {"run", "software_trigger_pre_ext", "software_trigger_post_ext"});
    } else if (node.HasColumn("software_trigger_pre")) {
        node = node.Define(
            "software_trigger",
            [](unsigned run, int pre, int post) {
                return run < 16880 ? pre > 0 : post > 0;
            },
            {"run", "software_trigger_pre", "software_trigger_post"});
    }

    return node;
}

const rarexsec::Processor& rarexsec::processor() {
    static const Processor ep{};
    return ep;
}
