#pragma once
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <TMatrixDSym.h>
#include <TH1D.h>

#include "rarexsec/Hub.hh"
#include "rarexsec/plot/Plotter.hh"
#include "rarexsec/syst/Systematics.hh"

namespace rarexsec::syst {

// Beam ⊕ Strangeness inputs and optional data
struct UnionSamples {
    std::vector<const Entry*> A_beam;        // kind == beam (non-strange)
    std::vector<const Entry*> B_strange;     // kind == strangeness
    std::vector<const Entry*> data;          // kind == data
};

// Physics/source configuration
struct UnionConfig {
    bool use_stat   = true;
    bool use_ppfx   = true;     // NuMI PPFX: weightsPPFX (+ ppfx_cv)
    bool use_genie  = true;     // map key "All_UBGenie" (+ weightSplineTimesTune)
    bool use_reint  = true;     // map key "reint_all"
    bool use_pot    = true;     // global norm, fully correlated A & B
    bool use_detvar = true;     // detector variations

    // Universe counts: -1 -> auto-detect from ntuples
    int N_ppfx   = -1;          // typically 600 (NuMI)
    int N_genie  = -1;          // you write 500
    int N_reint  = -1;          // often 100

    // Branch/key names (align with EventWeightAnalysis)
    std::string ppfx_branch     = "weightsPPFX";
    std::string ppfx_cv_branch  = "ppfx_cv";
    std::string map_branch      = "weights";
    std::string genie_key       = "All_UBGenie";
    std::string genie_cv_branch = "weightSplineTimesTune";
    std::string reint_key       = "reint_all";

    // Detector variations
    std::vector<std::pair<std::string,std::string>> detvar_pairs; // ± tags
    std::vector<std::string> detvar_unisims;                      // single-sided tags

    // POT fractional uncertainty (fully correlated across A & B)
    double pot_frac = 0.0;

    // Also build the A+B summed spectrum/covariances
    bool make_sum = true;
};

// Outputs: nominal hists and covariances (block and sum)
struct UnionProducts {
    std::unique_ptr<TH1D> H_A;
    std::unique_ptr<TH1D> H_B;
    std::unique_ptr<TH1D> H_sum;    // if make_sum
    std::unique_ptr<TH1D> H_data;   // optional

    TMatrixDSym C_block_total;
    std::map<std::string, TMatrixDSym> C_block_sources;

    TMatrixDSym C_sum_total;
    std::map<std::string, TMatrixDSym> C_sum_sources;
};

// Collect A, B, data from Hub
UnionSamples collect_union_samples(const Hub& hub,
                                   const std::string& beamline,
                                   const std::vector<std::string>& periods);

// Auto-detect universe counts (returns >0 if found, else default_val)
int detect_n_univ_ushort(const plot::H1Spec& spec,
                         const std::vector<const Entry*>& mc,
                         const std::string& branch,
                         int default_val = 0);

int detect_n_univ_map(const plot::H1Spec& spec,
                      const std::vector<const Entry*>& mc,
                      const std::string& map_branch,
                      const std::string& key,
                      int default_val = 0);

// Build covariances for a single spec used in both categories (same binning)
UnionProducts build_union_systematics(const plot::H1Spec& spec,
                                      const UnionSamples& samp,
                                      const UnionConfig& cfg);

// One-call runner that also collects samples from the Hub
UnionProducts run_union_systematics(const Hub& hub,
                                    const std::string& beamline,
                                    const std::vector<std::string>& periods,
                                    const plot::H1Spec& spec,
                                    const UnionConfig& cfg);

// Utility: union of A and B entries (for stacked plots)
std::vector<const Entry*> union_mc(const UnionSamples& s);

} // namespace rarexsec::syst
