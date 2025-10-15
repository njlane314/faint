#include "rarexsec/syst/SystematicsPack.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <cmath>
#include <stdexcept>

#include "rarexsec/proc/DataModel.h"
#include "rarexsec/syst/Systematics.h"

namespace rarexsec::systpack {

namespace {
std::unique_ptr<TH1D> clone_reset_like(const TH1D& templ, const std::string& name) {
  auto h = std::unique_ptr<TH1D>(static_cast<TH1D*>(templ.Clone(name.c_str())));
  h->SetDirectory(nullptr);
  h->Reset("ICES");
  return h;
}
//_______________________________________________________________________________________
std::unique_ptr<TH1D> sum_parts(std::vector<ROOT::RDF::RResultPtr<TH1D>>& parts,
                                const TH1D& model, const std::string& name) 
{
  std::unique_ptr<TH1D> total;
  for (auto& rr : parts) {
    const TH1D& h = rr.GetValue();
    if (!total) {
      total.reset(static_cast<TH1D*>(h.Clone(name.c_str())));
      total->SetDirectory(nullptr);
    } else {
      total->Add(&h);
    }
  }
  if (!total) return clone_reset_like(model, name);
  return total;
}
//_______________________________________________________________________________________
std::unique_ptr<TH1D> make_total_hist(const TH1D& model,
                                      const std::string& value_col,
                                      const std::string& weight_col,
                                      const std::vector<const rarexsec::Entry*>& entries,
                                      const std::string& name_suffix) 
{
  std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
  parts.reserve(entries.size());
  for (auto* e : entries) {
    if (!e) continue;
    auto node = e->rnode();
    parts.emplace_back(node.Histo1D(model, value_col, weight_col));
  }
  return sum_parts(parts, model, std::string(model.GetName()) + name_suffix);
}
//_______________________________________________________________________________________
std::unique_ptr<TH1D> make_total_hist_universe_ushort(const TH1D& model,
                                                      const std::string& value_col,
                                                      const std::string& base_weight_col,
                                                      const std::vector<const rarexsec::Entry*>& entries,
                                                      const std::string& weights_branch,
                                                      int k,
                                                      double us_scale,
                                                      const std::string& cv_branch,
                                                      const std::string& name_suffix) 
{
  std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
  parts.reserve(entries.size());
  for (size_t ie = 0; ie < entries.size(); ++ie) {
    auto* e = entries[ie];
    if (!e) continue;
    auto node = e->rnode();
    const std::string col = "_rx_univ_" + std::to_string(k) + "_src" + std::to_string(ie);
    if (cv_branch.empty()) {
      auto n1 = node.Define(col,
        [k, us_scale](const ROOT::RVec<unsigned short>& v, double w_nom) {
          double wk = 1.0;
          if (k >= 0 && k < (int)v.size()) wk = static_cast<double>(v[k]) * us_scale;
          const double out = w_nom * wk;
          return (std::isfinite(out) && out > 0.0) ? out : 0.0;
        }, {weights_branch, base_weight_col});
      parts.emplace_back(n1.Histo1D(model, value_col, col));
    } else {
      auto n1 = node.Define(col,
        [k, us_scale](const ROOT::RVec<unsigned short>& v, double w_nom, double w_cv) {
          double wk = 1.0;
          if (k >= 0 && k < (int)v.size()) wk = static_cast<double>(v[k]) * us_scale;
          const double out = w_nom * w_cv * wk;
          return (std::isfinite(out) && out > 0.0) ? out : 0.0;
        }, {weights_branch, base_weight_col, cv_branch});
      parts.emplace_back(n1.Histo1D(model, value_col, col));
    }
  }
  return sum_parts(parts, model, std::string(model.GetName()) + name_suffix);
}
} 
//_______________________________________________________________________________________
SystematicsPack::SystematicsPack(Config cfg) : cfg_(std::move(cfg)) {}
//_______________________________________________________________________________________
Result SystematicsPack::build(const TH1D& model,
                              const std::vector<const rarexsec::Entry*>& mc_entries,
                              const std::vector<const rarexsec::Entry*>& ext_entries) const 
{
  using rarexsec::syst::mc_stat_covariance;
  using rarexsec::syst::sample_covariance;
  using rarexsec::syst::sum;

  Result out;

  auto H_mc = make_total_hist(model, cfg_.value_col, cfg_.weight_col, mc_entries, "_mc");
  if (!H_mc) throw std::runtime_error("SystematicsPack: MC nominal is empty");

  out.sources["MC stat"] = mc_stat_covariance(*H_mc);

  if (cfg_.use_ppfx && cfg_.N_ppfx > 0) {
    std::vector<std::unique_ptr<TH1D>> universes;
    universes.reserve(cfg_.N_ppfx);
    for (int k = 0; k < cfg_.N_ppfx; ++k) {
      universes.emplace_back(
        make_total_hist_universe_ushort(model, cfg_.value_col, cfg_.weight_col,
                                        mc_entries, cfg_.ppfx_branch, k,
                                        cfg_.ushort_scale, cfg_.ppfx_cv_branch,
                                        "_ppfx_u" + std::to_string(k)));
    }
    out.sources["Flux (PPFX)"] = sample_covariance(*H_mc, universes);
  }

  if (cfg_.use_genie && cfg_.N_genie > 0) {
    std::vector<std::unique_ptr<TH1D>> universes;
    universes.reserve(cfg_.N_genie);
    for (int k = 0; k < cfg_.N_genie; ++k) {
      universes.emplace_back(
        make_total_hist_universe_ushort(model, cfg_.value_col, cfg_.weight_col,
                                        mc_entries, cfg_.genie_branch, k,
                                        cfg_.ushort_scale, cfg_.genie_cv_branch,
                                        "_genie_u" + std::to_string(k)));
    }
    out.sources["GENIE"] = sample_covariance(*H_mc, universes);
  }

  if (cfg_.use_reint && cfg_.N_reint > 0) {
    std::vector<std::unique_ptr<TH1D>> universes;
    universes.reserve(cfg_.N_reint);
    for (int k = 0; k < cfg_.N_reint; ++k) {
      universes.emplace_back(
        make_total_hist_universe_ushort(model, cfg_.value_col, cfg_.weight_col,
                                        mc_entries, cfg_.reint_branch, k,
                                        cfg_.ushort_scale, "",
                                        "_reint_u" + std::to_string(k)));
    }
    out.sources["Reint (Geant4)"] = sample_covariance(*H_mc, universes);
  }

  out.H_pred = std::unique_ptr<TH1D>(static_cast<TH1D*>(H_mc->Clone("H_pred")));
  out.H_pred->SetDirectory(nullptr);

  if (cfg_.include_ext && !ext_entries.empty()) {
    if (auto H_ext = make_total_hist(model, cfg_.value_col, cfg_.weight_col, ext_entries, "_ext")) {
      out.H_pred->Add(H_ext.get());
      out.sources["EXT stat"] = mc_stat_covariance(*H_ext);
    }
  }

  std::vector<const TMatrixDSym*> pieces;
  pieces.reserve(out.sources.size());
  for (auto& kv : out.sources) pieces.push_back(&kv.second);
  out.total = sum(pieces);

  return out;
}
//_______________________________________________________________________________________
} // namespace rarexsec::systpack
