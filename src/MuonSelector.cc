#include "faint/MuonSelector.h"

#include <cmath>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

#include "ROOT/RVec.hxx"
#include "faint/FiducialVolume.h"
#include "faint/Selections.h"

namespace faint {

ROOT::RDF::RNode MuonSelector::process(ROOT::RDF::RNode df,
                                       SampleOrigin origin) const {
  if (!df.HasColumn("track_shower_scores")) {
    auto no_mu_df = df.Define("n_muons_tot", []() { return 0UL; })
                        .Define("has_muon", []() { return false; });
    return next_ ? next_->process(no_mu_df, origin) : no_mu_df;
  }

  auto muon_mask_df = this->build_mask(df);
  auto muon_features_df = this->extract_features(muon_mask_df);
  return next_ ? next_->process(muon_features_df, origin) : muon_features_df;
}

ROOT::RDF::RNode MuonSelector::build_mask(ROOT::RDF::RNode df) const {
  return df.Define(
      "muon_mask",
      [](const ROOT::RVec<float> &scores, const ROOT::RVec<float> &llr,
         const ROOT::RVec<float> &lengths, const ROOT::RVec<float> &dists,
         const ROOT::RVec<float> &start_x, const ROOT::RVec<float> &start_y,
         const ROOT::RVec<float> &start_z, const ROOT::RVec<float> &end_x,
         const ROOT::RVec<float> &end_y, const ROOT::RVec<float> &end_z,
         const ROOT::RVec<unsigned> &gens) {
        ROOT::RVec<bool> mask(scores.size());
        for (size_t i = 0; i < scores.size(); ++i) {
          const bool fid_start = fiducial::is_in_reco_volume(
              start_x[i], start_y[i], start_z[i]);
          const bool fid_end = fiducial::is_in_reco_volume(
              end_x[i], end_y[i], end_z[i]);
          mask[i] = selection::passes_muon_track_selection(
              scores[i], llr[i], lengths[i], dists[i], gens[i], fid_start,
              fid_end);
        }
        return mask;
      },
      {"track_shower_scores", "trk_llr_pid_v", "track_length",
       "track_distance_to_vertex", "track_start_x", "track_start_y",
       "track_start_z", "track_end_x", "track_end_y", "track_end_z",
       "pfp_generations"});
}

ROOT::RDF::RNode MuonSelector::extract_features(ROOT::RDF::RNode df) const {
  auto filter_float = [](const ROOT::RVec<float> &vals,
                         const ROOT::RVec<bool> &mask) {
    ROOT::RVec<float> out;
    out.reserve(vals.size());
    for (size_t i = 0; i < vals.size(); ++i)
      if (mask[i])
        out.push_back(vals[i]);
    return out;
  };

  auto filter_uint = [](const ROOT::RVec<unsigned> &vals,
                        const ROOT::RVec<bool> &mask) {
    ROOT::RVec<unsigned> out;
    out.reserve(vals.size());
    for (size_t i = 0; i < vals.size(); ++i)
      if (mask[i])
        out.push_back(vals[i]);
    return out;
  };

  auto filter_costheta = [](const ROOT::RVec<float> &theta,
                            const ROOT::RVec<bool> &mask) {
    ROOT::RVec<float> out;
    out.reserve(theta.size());
    for (size_t i = 0; i < theta.size(); ++i)
      if (mask[i])
        out.push_back(std::cos(theta[i]));
    return out;
  };

  auto mu_df = df.Define("muon_trk_score_v", filter_float,
                         {"track_shower_scores", "muon_mask"})
                   .Define("muon_trk_llr_pid_v", filter_float,
                           {"trk_llr_pid_v", "muon_mask"})
                   .Define("muon_trk_start_x_v", filter_float,
                           {"track_start_x", "muon_mask"})
                   .Define("muon_trk_start_y_v", filter_float,
                           {"track_start_y", "muon_mask"})
                   .Define("muon_trk_start_z_v", filter_float,
                           {"track_start_z", "muon_mask"})
                   .Define("muon_trk_end_x_v", filter_float,
                           {"track_end_x", "muon_mask"})
                   .Define("muon_trk_end_y_v", filter_float,
                           {"track_end_y", "muon_mask"})
                   .Define("muon_trk_end_z_v", filter_float,
                           {"track_end_z", "muon_mask"})
                   .Define("muon_trk_length_v", filter_float,
                           {"track_length", "muon_mask"})
                   .Define("muon_trk_distance_v", filter_float,
                           {"track_distance_to_vertex", "muon_mask"})
                   .Define("muon_pfp_generation_v", filter_uint,
                           {"pfp_generations", "muon_mask"})
                   .Define("muon_track_costheta", filter_costheta,
                           {"track_theta", "muon_mask"});

  auto redefine_or_define = [](ROOT::RDF::RNode node, const char *name,
                               auto &&callable,
                               std::initializer_list<std::string> cols) {
    const std::vector<std::string> columns(cols.begin(), cols.end());
    return node.HasColumn(name) ? node.Redefine(name, std::forward<decltype(callable)>(callable), columns)
                                : node.Define(name, std::forward<decltype(callable)>(callable), columns);
  };

  mu_df = redefine_or_define(
      mu_df, "n_muons_tot",
      [](const ROOT::RVec<bool> &mask) {
        return static_cast<unsigned long>(ROOT::VecOps::Sum(mask));
      },
      {"muon_mask"});

  mu_df = redefine_or_define(
      mu_df, selection::column::kPassMuon,
      [](unsigned long n_muons) {
        return selection::passes_muon_selection(n_muons);
      },
      {"n_muons_tot"});

  mu_df = redefine_or_define(
      mu_df, selection::column::kPassFinal,
      [](bool pre, bool flash, bool fiducial, bool muon, bool topo) {
        return selection::passes_final_selection(pre, flash, fiducial, muon,
                                                 topo);
      },
      {selection::column::kPassPre, selection::column::kPassFlash,
       selection::column::kPassFiducial, selection::column::kPassMuon,
       selection::column::kPassTopology});

  mu_df = redefine_or_define(
      mu_df, "has_muon",
      [](unsigned long n_muons) {
        return selection::passes_muon_selection(n_muons);
      },
      {"n_muons_tot"});
  return mu_df;
}

} // namespace faint
