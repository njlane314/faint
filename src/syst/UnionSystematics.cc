#include "rarexsec/syst/UnionSystematics.hh"

#include <ROOT/RDataFrame.hxx>
#include <stdexcept>

namespace rarexsec::syst {

using MapSD = std::map<std::string, std::vector<double>>;

UnionSamples collect_union_samples(const Hub& hub,
                                   const std::string& beamline,
                                   const std::vector<std::string>& periods) {
    UnionSamples s;
    const auto sim  = hub.simulation_entries(beamline, periods);
    const auto data = hub.data_entries(beamline, periods);

    for (auto* e : sim) {
        if (!e) continue;
        if (e->kind == sample::origin::beam)            s.A_beam.push_back(e);
        else if (e->kind == sample::origin::strangeness) s.B_strange.push_back(e);
    }
    s.data = data;
    return s;
}

int detect_n_univ_ushort(const plot::H1Spec& spec,
                         const std::vector<const Entry*>& mc,
                         const std::string& branch,
                         int default_val) {
    for (auto* e : mc) {
        if (!e) continue;
        auto n0 = selection::apply(e->rnode(), spec.sel, *e);
        auto n1 = n0.Define("_rx_len_", [](const ROOT::RVec<unsigned short>& v){ return int(v.size()); },
                            {branch}).Range(1);
        auto lens = n1.Take<int>("_rx_len_");
        if (!lens->empty() && (*lens)[0] > 0) return (*lens)[0];
    }
    return default_val;
}

int detect_n_univ_map(const plot::H1Spec& spec,
                      const std::vector<const Entry*>& mc,
                      const std::string& map_branch,
                      const std::string& key,
                      int default_val) {
    for (auto* e : mc) {
        if (!e) continue;
        auto n0 = selection::apply(e->rnode(), spec.sel, *e);
        auto n1 = n0.Define("_rx_len_", [key](const MapSD& m){
                                auto it = m.find(key);
                                return (it != m.end()) ? int(it->second.size()) : 0;
                            }, {map_branch}).Range(1);
        auto lens = n1.Take<int>("_rx_len_");
        if (!lens->empty() && (*lens)[0] > 0) return (*lens)[0];
    }
    return default_val;
}

std::vector<const Entry*> union_mc(const UnionSamples& s) {
    std::vector<const Entry*> out; out.reserve(s.A_beam.size()+s.B_strange.size());
    out.insert(out.end(), s.A_beam.begin(), s.A_beam.end());
    out.insert(out.end(), s.B_strange.begin(), s.B_strange.end());
    return out;
}

// Build block-level covariance from single-sided detvar tags by sample covariance over tags
static TMatrixDSym block_cov_from_detvar_unisims(const plot::H1Spec& specA, const std::vector<const Entry*>& A,
                                                 const plot::H1Spec& specB, const std::vector<const Entry*>& B,
                                                 const std::vector<std::string>& tags) {
    if (tags.empty()) return TMatrixDSym(0);
    auto H0A = make_total_mc_hist(specA, A, "_A_nom");
    auto H0B = make_total_mc_hist(specB, B, "_B_nom");
    if (!H0A || !H0B) return TMatrixDSym(0);

    const int nA = H0A->GetNbinsX();
    const int nB = H0B->GetNbinsX();
    const int n  = nA + nB;

    std::vector<std::vector<double>> deltas; deltas.reserve(tags.size());
    for (const auto& t : tags) {
        auto HA = make_total_mc_hist_detvar(specA, A, t, "_A");
        auto HB = make_total_mc_hist_detvar(specB, B, t, "_B");
        if (!HA || !HB) continue;

        std::vector<double> d(n, 0.0);
        for (int i=0;i<nA;++i) d[i]      = HA->GetBinContent(i+1) - H0A->GetBinContent(i+1);
        for (int j=0;j<nB;++j) d[nA + j] = HB->GetBinContent(j+1) - H0B->GetBinContent(j+1);
        deltas.emplace_back(std::move(d));
    }

    const int N = (int)deltas.size();
    if (N <= 1) return TMatrixDSym(n);

    TMatrixDSym C(n);
    for (int i=0;i<n;++i) {
        for (int j=i;j<n;++j) {
            long double s = 0.0L;
            for (int u=0; u<N; ++u) s += (long double)deltas[u][i] * (long double)deltas[u][j];
            double cij = (double)(s / (N - 1));
            C(i,j) = C(j,i) = cij;
        }
    }
    return C;
}

UnionProducts build_union_systematics(const plot::H1Spec& spec,
                                      const UnionSamples& samp,
                                      const UnionConfig& cfg) {
    UnionProducts out;

    // Nominals
    out.H_A = make_total_mc_hist(spec, samp.A_beam, "_A");
    out.H_B = make_total_mc_hist(spec, samp.B_strange, "_B");
    if (!out.H_A || !out.H_B) throw std::runtime_error("UnionSystematics: empty A or B nominal");
    if (!samp.data.empty()) out.H_data = make_total_mc_hist(spec, samp.data, "_data");

    const int nA = out.H_A->GetNbinsX();
    const int nB = out.H_B->GetNbinsX();

    // Sources (block level on [A;B])
    if (cfg.use_stat) {
        out.C_block_sources["MC stat"] = block_diag_stat(*out.H_A, *out.H_B);
    }

    if (cfg.use_ppfx) {
        int N = cfg.N_ppfx;
        if (N < 0) N = detect_n_univ_ushort(spec, samp.A_beam, cfg.ppfx_branch, 600);
        out.C_block_sources["Flux (PPFX)"] =
            block_cov_from_weight_vector_ushort_scaled(
                spec, samp.A_beam, spec, samp.B_strange,
                cfg.ppfx_branch, N, cfg.ppfx_cv_branch, 1.0/1000.0);
    }

    if (cfg.use_genie) {
        int N = cfg.N_genie;
        if (N < 0) N = detect_n_univ_map(spec, samp.A_beam, cfg.map_branch, cfg.genie_key, 500);
        out.C_block_sources["GENIE"] =
            block_cov_from_map_weight_vector(
                spec, samp.A_beam, spec, samp.B_strange,
                cfg.map_branch, cfg.genie_key, N, cfg.genie_cv_branch);
    }

    if (cfg.use_reint) {
        int N = cfg.N_reint;
        if (N < 0) N = detect_n_univ_map(spec, samp.A_beam, cfg.map_branch, cfg.reint_key, 100);
        out.C_block_sources["Reint (Geant4)"] =
            block_cov_from_map_weight_vector(
                spec, samp.A_beam, spec, samp.B_strange,
                cfg.map_branch, cfg.reint_key, N, "");
    }

    if (cfg.use_pot && cfg.pot_frac > 0.0) {
        out.C_block_sources["POT"] = pot_cov_block(*out.H_A, *out.H_B, cfg.pot_frac);
    }

    if (cfg.use_detvar) {
        TMatrixDSym Cdet(out.H_A->GetNbinsX() + out.H_B->GetNbinsX());
        bool any = false;

        if (!cfg.detvar_pairs.empty()) {
            Cdet += block_cov_from_detvar_pairs(
                spec, samp.A_beam, spec, samp.B_strange, cfg.detvar_pairs);
            any = true;
        }
        if (!cfg.detvar_unisims.empty()) {
            Cdet += block_cov_from_detvar_unisims(
                spec, samp.A_beam, spec, samp.B_strange, cfg.detvar_unisims);
            any = true;
        }
        if (any) out.C_block_sources["Detector"] = std::move(Cdet);
    }

    // Sum block → total block
    {
        std::vector<const TMatrixDSym*> pieces;
        pieces.reserve(out.C_block_sources.size());
        for (auto& kv : out.C_block_sources) pieces.push_back(&kv.second);
        out.C_block_total = sum(pieces);
    }

    // Map to A+B summed spectrum (same binning spec for A and B)
    if (cfg.make_sum) {
        out.H_sum = sum_same_binning(*out.H_A, *out.H_B, "h_sum");
        for (auto& kv : out.C_block_sources) {
            out.C_sum_sources[kv.first] =
                sum_covariance_block_same_binning(kv.second, nA, nB);
        }
        out.C_sum_total = sum_covariance_block_same_binning(out.C_block_total, nA, nB);
    }

    return out;
}

UnionProducts run_union_systematics(const Hub& hub,
                                    const std::string& beamline,
                                    const std::vector<std::string>& periods,
                                    const plot::H1Spec& spec,
                                    const UnionConfig& cfg) {
    auto samples = collect_union_samples(hub, beamline, periods);
    return build_union_systematics(spec, samples, cfg);
}

} // namespace rarexsec::syst
