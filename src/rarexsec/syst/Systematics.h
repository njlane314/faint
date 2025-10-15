#pragma once
#include <TH1D.h>
#include <TMatrixDSym.h>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "rarexsec/Hub.h"
#include "rarexsec/plot/Plotter.h"
#include "rarexsec/proc/Selection.h"

namespace rarexsec::syst {

inline constexpr int RAREXSEC_MULTISIM_DDOF = 1;

TMatrixDSym mc_stat_covariance(const TH1D&);
TMatrixDSym sample_covariance(const TH1D&, const std::vector<std::unique_ptr<TH1D>>&);
TMatrixDSym hessian_covariance(const TH1D&, const TH1D&, const TH1D&);
TMatrixDSym sum(const std::vector<const TMatrixDSym*>&);
std::unique_ptr<TH1D> make_total_mc_hist(const plot::TH1DModel& spec,
                                         const std::vector<const Entry*>& entries,
                                         const std::string& suffix);
std::unique_ptr<TH1D> make_total_mc_hist_detvar(const plot::TH1DModel& spec,
                                                const std::vector<const Entry*>& entries,
                                                const std::string& tag,
                                                const std::string& suffix);

TMatrixDSym cov_from_detvar_pairs(
    const plot::TH1DModel& spec, const std::vector<const Entry*>& mc,
    const std::vector<std::pair<std::string, std::string>>& tag_pairs);

TMatrixDSym cov_from_detvar_unisims(
    const plot::TH1DModel& spec, const std::vector<const Entry*>& mc,
    const std::vector<std::string>& tags);

TMatrixDSym block_cov_from_detvar_pairs(
    const plot::TH1DModel& specA, const std::vector<const Entry*>& A,
    const plot::TH1DModel& specB, const std::vector<const Entry*>& B,
    const std::vector<std::pair<std::string, std::string>>& tag_pairs);

std::unique_ptr<TH1D> make_total_mc_hist_weight_universe_ushort(
    const plot::TH1DModel& spec, const std::vector<const Entry*>& mc,
    const std::string& weights_branch, int k, const std::string& suffix,
    const std::string& cv_branch = "", double us_scale = 1.0 / 1000.0);

TMatrixDSym cov_from_weight_vector_ushort(
    const plot::TH1DModel& spec, const std::vector<const Entry*>& mc,
    const std::string& weights_branch, int nuniv,
    const std::string& cv_branch = "", double us_scale = 1.0 / 1000.0);

std::unique_ptr<TH1D> make_total_mc_hist_weight_universe_map(
    const plot::TH1DModel& spec, const std::vector<const Entry*>& mc,
    const std::string& map_branch, const std::string& key, int k,
    const std::string& suffix, const std::string& cv_branch = "");

TMatrixDSym cov_from_map_weight_vector(
    const plot::TH1DModel& spec, const std::vector<const Entry*>& mc,
    const std::string& map_branch, const std::string& key, int nuniv,
    const std::string& cv_branch = "");

TMatrixDSym block_cov_from_weight_vector_ushort_scaled(
    const plot::TH1DModel& specA, const std::vector<const Entry*>& A,
    const plot::TH1DModel& specB, const std::vector<const Entry*>& B,
    const std::string& weights_branch, int nuniv,
    const std::string& cv_branch = "", double us_scale = 1.0 / 1000.0);

TMatrixDSym block_cov_from_map_weight_vector(
    const plot::TH1DModel& specA, const std::vector<const Entry*>& A,
    const plot::TH1DModel& specB, const std::vector<const Entry*>& B,
    const std::string& map_branch, const std::string& key, int nuniv,
    const std::string& cv_branch = "");

TMatrixDSym block_cov_from_ud_ushort(
    const plot::TH1DModel& specA, const std::vector<const Entry*>& A,
    const plot::TH1DModel& specB, const std::vector<const Entry*>& B,
    const std::string& up_branch, const std::string& dn_branch, int knob_index,
    double us_scale = 1.0 / 1000.0, const std::string& cv_branch = "");

TMatrixDSym block_diag_stat(const TH1D& A, const TH1D& B);
TMatrixDSym pot_cov_block(const TH1D& A, const TH1D& B, double frac_pot);

std::unique_ptr<TH1D> sum_same_binning(const TH1D& A, const TH1D& B, const std::string& name);

TMatrixDSym sum_covariance_block_same_binning(const TMatrixDSym& C_block, int nA, int nB);

}
