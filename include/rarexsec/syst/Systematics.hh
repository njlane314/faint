#pragma once
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <TMatrixDSym.h>
#include <TH1D.h>

#include "rarexsec/Plotter.hh"
#include "rarexsec/Hub.hh"
#include "rarexsec/proc/Selection.hh"

namespace rarexsec::syst {

// ---- Already existing APIs (yours) ----
// TMatrixDSym mc_stat_covariance(const TH1D&);
// TMatrixDSym sample_covariance(const TH1D&, const std::vector<std::unique_ptr<TH1D>>&);
// TMatrixDSym hessian_covariance(const TH1D&, const TH1D&, const TH1D&);
// TMatrixDSym sum(const std::vector<const TMatrixDSym*>&);
// TMatrixDSym shape_only(const TMatrixDSym&, const TH1D&);
// std::unique_ptr<TH1D> make_total_mc_hist(...);
// std::unique_ptr<TH1D> make_total_mc_hist_weight_universe(...);
// std::unique_ptr<TH1D> make_total_mc_hist_detvar(...);
// TMatrixDSym cov_from_weight_vector(...);
// TMatrixDSym cov_from_detvar_pm(...);
TMatrixDSym cov_from_detvar_pairs(
    const plot::H1Spec& spec, const std::vector<const Entry*>& mc,
    const std::vector<std::pair<std::string,std::string>>& tag_pairs);

TMatrixDSym cov_from_detvar_unisims(
    const plot::H1Spec& spec, const std::vector<const Entry*>& mc,
    const std::vector<std::string>& tags);

TMatrixDSym block_cov_from_detvar_pairs(
    const plot::H1Spec& specA, const std::vector<const Entry*>& A,
    const plot::H1Spec& specB, const std::vector<const Entry*>& B,
    const std::vector<std::pair<std::string,std::string>>& tag_pairs);

// ---- New: unsigned-short universe vectors (+ optional CV multiplier) ----
std::unique_ptr<TH1D> make_total_mc_hist_weight_universe_ushort(
    const plot::H1Spec& spec, const std::vector<const Entry*>& mc,
    const std::string& weights_branch, int k, const std::string& suffix,
    const std::string& cv_branch = "", double us_scale = 1.0/1000.0);

TMatrixDSym cov_from_weight_vector_ushort(
    const plot::H1Spec& spec, const std::vector<const Entry*>& mc,
    const std::string& weights_branch, int nuniv,
    const std::string& cv_branch = "", double us_scale = 1.0/1000.0);

// ---- New: map<string, vector<double>> universes ----
std::unique_ptr<TH1D> make_total_mc_hist_weight_universe_map(
    const plot::H1Spec& spec, const std::vector<const Entry*>& mc,
    const std::string& map_branch, const std::string& key, int k,
    const std::string& suffix, const std::string& cv_branch = "");

TMatrixDSym cov_from_map_weight_vector(
    const plot::H1Spec& spec, const std::vector<const Entry*>& mc,
    const std::string& map_branch, const std::string& key, int nuniv,
    const std::string& cv_branch = "");

// ---- New: block covariances for (A=beam ⊕ B=strangeness) ----
TMatrixDSym block_cov_from_weight_vector_ushort_scaled(
    const plot::H1Spec& specA, const std::vector<const Entry*>& A,
    const plot::H1Spec& specB, const std::vector<const Entry*>& B,
    const std::string& weights_branch, int nuniv,
    const std::string& cv_branch = "", double us_scale = 1.0/1000.0);

TMatrixDSym block_cov_from_map_weight_vector(
    const plot::H1Spec& specA, const std::vector<const Entry*>& A,
    const plot::H1Spec& specB, const std::vector<const Entry*>& B,
    const std::string& map_branch, const std::string& key, int nuniv,
    const std::string& cv_branch = "");

TMatrixDSym block_cov_from_ud_ushort(
    const plot::H1Spec& specA, const std::vector<const Entry*>& A,
    const plot::H1Spec& specB, const std::vector<const Entry*>& B,
    const std::string& up_branch, const std::string& dn_branch, int knob_index,
    double us_scale = 1.0/1000.0, const std::string& cv_branch = "");

// ---- New: POT and block-diag stat helpers ----
TMatrixDSym block_diag_stat(const TH1D& A, const TH1D& B);
TMatrixDSym pot_cov_block(const TH1D& A, const TH1D& B, double frac_pot);

// ---- New: combine categories with same binning (sum A⊕B per bin) ----
std::unique_ptr<TH1D> sum_same_binning(const TH1D& A, const TH1D& B, const std::string& name);

// Map block covariance C_{(A⊕B)} to the summed spectrum covariance.
// nA = #bins of A; nB = #bins of B (must match A for this routine).
TMatrixDSym sum_covariance_block_same_binning(const TMatrixDSym& C_block, int nA, int nB);

} // namespace rarexsec::syst
