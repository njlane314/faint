#pragma once
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <ROOT/RDataFrame.hxx>
#include <TMatrixDSym.h>
#include <TH1D.h>

#include "rarexsec/Plotter.hh"
#include "rarexsec/Hub.hh"
#include "rarexsec/proc/Selection.hh"

namespace rarexsec::syst {

TMatrixDSym mc_stat_covariance(const TH1D& h);

TMatrixDSym sample_covariance(const TH1D& nominal,
                              const std::vector<std::unique_ptr<TH1D>>& universes);

TMatrixDSym hessian_covariance(const TH1D& nominal,
                               const TH1D& plus,
                               const TH1D& minus);

TMatrixDSym sum(const std::vector<const TMatrixDSym*>& terms);

TMatrixDSym shape_only(const TMatrixDSym& cov, const TH1D& nominal);

std::unique_ptr<TH1D> make_total_mc_hist(const plot::H1Spec& spec,
                                         const std::vector<const Entry*>& mc,
                                         const std::string& suffix = "_nominal");

std::unique_ptr<TH1D> make_total_mc_hist_weight_universe(const plot::H1Spec& spec,
                                                         const std::vector<const Entry*>& mc,
                                                         const std::string& weights_branch,
                                                         int k,
                                                         const std::string& suffix);

std::unique_ptr<TH1D> make_total_mc_hist_detvar(const plot::H1Spec& spec,
                                                const std::vector<const Entry*>& mc,
                                                const std::string& detvar_tag,
                                                const std::string& suffix);

TMatrixDSym cov_from_weight_vector(const plot::H1Spec& spec,
                                   const std::vector<const Entry*>& mc,
                                   const std::string& weights_branch,
                                   int nuniv);

TMatrixDSym cov_from_detvar_pm(const plot::H1Spec& spec,
                               const std::vector<const Entry*>& mc,
                               const std::string& tag_up,
                               const std::string& tag_down);

}
