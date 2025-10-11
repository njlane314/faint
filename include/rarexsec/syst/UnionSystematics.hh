#pragma once
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <TMatrixDSym.h>
#include <TH1D.h>

#include "rarexsec/Hub.hh"
#include "rarexsec/Plotter.hh"
#include "rarexsec/syst/Systematics.hh"

namespace rarexsec::syst {

struct UnionSamples {
    std::vector<const Entry*> A_beam;
    std::vector<const Entry*> B_strange;
    std::vector<const Entry*> dirt;
    std::vector<const Entry*> ext;
    std::vector<const Entry*> data;
};

struct UnionConfig {
    bool use_stat   = true;
    bool use_ppfx   = true;
    bool use_genie  = true;
    bool use_reint  = true;
    bool use_pot    = true;
    bool use_detvar = true;

    bool include_dirt = true;
    bool include_ext  = true;

    int N_ppfx  = -1;
    int N_genie = -1;
    int N_reint = -1;

    std::string ppfx_branch     = "weightsPPFX";
    std::string ppfx_cv_branch  = "ppfx_cv";
    std::string map_branch      = "weights";
    std::string genie_key       = "All_UBGenie";
    std::string genie_cv_branch = "weightSplineTimesTune";
    std::string reint_key       = "reint_all";

    std::vector<std::pair<std::string,std::string>> detvar_pairs;
    std::vector<std::string> detvar_unisims;

    double pot_frac       = 0.0;
    double dirt_norm_frac = 0.0;
    double ext_norm_frac  = 0.0;

    bool make_sum = true;
};

struct UnionProducts {
    std::unique_ptr<TH1D> H_A;
    std::unique_ptr<TH1D> H_B;
    std::unique_ptr<TH1D> H_DIRT;
    std::unique_ptr<TH1D> H_EXT;
    std::unique_ptr<TH1D> H_sum;
    std::unique_ptr<TH1D> H_data;

    TMatrixDSym C_block_total;
    std::map<std::string, TMatrixDSym> C_block_sources;

    TMatrixDSym C_sum_total;
    std::map<std::string, TMatrixDSym> C_sum_sources;
};

UnionSamples collect_union_samples(const Hub& hub,
                                   const std::string& beamline,
                                   const std::vector<std::string>& periods);

int detect_n_univ_ushort(const plot::Histogram1DSpec& spec,
                         const std::vector<const Entry*>& mc,
                         const std::string& branch,
                         int default_val = 0);

int detect_n_univ_map(const plot::Histogram1DSpec& spec,
                      const std::vector<const Entry*>& mc,
                      const std::string& map_branch,
                      const std::string& key,
                      int default_val = 0);

UnionProducts build_union_systematics(const plot::Histogram1DSpec& spec,
                                      const UnionSamples& samp,
                                      const UnionConfig& cfg);

UnionProducts run_union_systematics(const Hub& hub,
                                    const std::string& beamline,
                                    const std::vector<std::string>& periods,
                                    const plot::Histogram1DSpec& spec,
                                    const UnionConfig& cfg);

std::vector<const Entry*> mc_union_AB(const UnionSamples& s);
std::vector<const Entry*> mc_union_all(const UnionSamples& s);

}
