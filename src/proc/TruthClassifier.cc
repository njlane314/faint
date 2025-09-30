#include "faint/proc/TruthClassifier.h"

#include "faint/FiducialVolume.h"

#include <cmath>
#include <iostream>
#include <map>
#include <mutex>

namespace faint {

using sample::SampleOrigin;

ROOT::RDF::RNode TruthClassifier::process(ROOT::RDF::RNode df,
                                          SampleOrigin origin) const {
  if (origin != SampleOrigin::kMonteCarlo) {
    return this->process_non_mc(df, origin);
  }

  auto count_frame = this->define_counts(df);
  auto inclusive_frame = this->assign_inclusive_channels(count_frame);
  auto exclusive_frame = this->assign_exclusive_channels(inclusive_frame);
  auto channel_frame = this->assign_channel_definitions(exclusive_frame);

  return next_ ? next_->process(channel_frame, origin) : channel_frame;
}

ROOT::RDF::RNode TruthClassifier::process_non_mc(ROOT::RDF::RNode df,
                                                 SampleOrigin origin) const {
  auto mode_frame = df.Define("genie_int_mode", []() { return -1; });

  auto inclusive_frame = mode_frame.Define("incl_channel", [c = origin]() {
    if (c == SampleOrigin::kData)
      return 0;
    if (c == SampleOrigin::kExternal)
      return 1;
    if (c == SampleOrigin::kDirt)
      return 2;
    return 99;
  });

  auto inclusive_alias_frame =
      inclusive_frame.Define("inclusive_strange_channels", "incl_channel");

  auto exclusive_frame = inclusive_alias_frame.Define("excl_channel", [c = origin]() {
    if (c == SampleOrigin::kData)
      return 0;
    if (c == SampleOrigin::kExternal)
      return 1;
    if (c == SampleOrigin::kDirt)
      return 2;
    return 99;
  });

  auto exclusive_alias_frame =
      exclusive_frame.Define("exclusive_strange_channels", "excl_channel");

  auto channel_frame = exclusive_alias_frame.Define("channel_def", [c = origin]() {
    if (c == SampleOrigin::kData)
      return 0;
    if (c == SampleOrigin::kExternal || c == SampleOrigin::kDirt)
      return 1;
    return 99;
  });

  auto channel_alias_frame = channel_frame.Define("channel_definitions", "channel_def");

  return next_ ? next_->process(channel_alias_frame, origin) : channel_alias_frame;
}

ROOT::RDF::RNode TruthClassifier::define_counts(ROOT::RDF::RNode df) const {
  auto fiducial_frame = df.Define(
      "in_fiducial",
      [](float x, float y, float z) {
        return fiducial::is_in_truth_volume(x, y, z);
      },
      {"neutrino_vertex_x", "neutrino_vertex_y", "neutrino_vertex_z"});

  auto strange_frame = fiducial_frame.Define(
      "mc_n_strange",
      "count_kaon_plus + count_kaon_minus + count_kaon_zero +"
      " count_lambda + count_sigma_plus + count_sigma_zero + count_sigma_minus");

  auto pion_frame = strange_frame.Define("mc_n_pion",
                                   "count_pi_plus + count_pi_minus");
  auto proton_frame = pion_frame.Define("mc_n_proton", "count_proton");

  auto mode_frame = proton_frame.Define(
      "genie_int_mode",
      [](int mode) {
        struct ModeCounter {
          std::map<int, long long> counts;
          std::mutex mtx;
          ~ModeCounter() {
            std::cout << "[DEBUG] GENIE interaction mode frequencies:\n";
            for (const auto &kv : counts) {
              std::cout << "  mode " << kv.first << ": " << kv.second
                        << std::endl;
            }
          }
        };
        static ModeCounter counter;
        {
          std::lock_guard<std::mutex> lock(counter.mtx);
          counter.counts[mode]++;
          if (counter.counts[mode] == 1 && mode != 0 && mode != 1 &&
              mode != 2 && mode != 3 && mode != 10) {
            std::cout << "[DEBUG] Uncategorised GENIE mode: " << mode
                      << std::endl;
          }
        }
        switch (mode) {
        case 0:
          return 0;
        case 1:
          return 1;
        case 2:
          return 2;
        case 3:
          return 3;
        case 10:
          return 10;
        default:
          return -1;
        }
      },
      {"interaction_mode"});

  return mode_frame;
}

ROOT::RDF::RNode
TruthClassifier::assign_inclusive_channels(ROOT::RDF::RNode df) const {
  auto inclusive_channel_frame = df.Define(
      "incl_channel",
      [](bool fv, int nu, int ccnc, int s, int npi, int np) {
        if (!fv)
          return 98;
        if (ccnc == 1)
          return 31;
        if (std::abs(nu) == 12 && ccnc == 0)
          return 30;
        if (std::abs(nu) == 14 && ccnc == 0) {
          if (s == 1)
            return 10;
          if (s > 1)
            return 11;
          if (np >= 1 && npi == 0)
            return 20;
          if (np == 0 && npi >= 1)
            return 21;
          if (np >= 1 && npi >= 1)
            return 22;
          return 23;
        }
        return 99;
      },
      {"in_fiducial", "neutrino_pdg", "interaction_ccnc", "mc_n_strange",
       "mc_n_pion", "mc_n_proton"});

  auto inclusive_alias_frame =
      inclusive_channel_frame.Define("inclusive_strange_channels", "incl_channel");
  return inclusive_alias_frame;
}

ROOT::RDF::RNode
TruthClassifier::assign_exclusive_channels(ROOT::RDF::RNode df) const {
  auto exclusive_channel_frame = df.Define(
      "excl_channel",
      [](bool fv, int nu, int ccnc, int s, int kp, int km, int k0, int lam,
         int sp, int s0, int sm) {
        if (!fv)
          return 98;
        if (ccnc == 1)
          return 31;
        if (std::abs(nu) == 12 && ccnc == 0)
          return 30;
        if (std::abs(nu) == 14 && ccnc == 0) {
          if (s == 0)
            return 32;
          if ((kp == 1 || km == 1) && s == 1)
            return 50;
          if (k0 == 1 && s == 1)
            return 51;
          if (lam == 1 && s == 1)
            return 52;
          if ((sp == 1 || sm == 1) && s == 1)
            return 53;
          if (lam == 1 && (kp == 1 || km == 1) && s == 2)
            return 54;
          if ((sp == 1 || sm == 1) && k0 == 1 && s == 2)
            return 55;
          if ((sp == 1 || sm == 1) && (kp == 1 || km == 1) && s == 2)
            return 56;
          if (lam == 1 && k0 == 1 && s == 2)
            return 57;
          if (kp == 1 && km == 1 && s == 2)
            return 58;
          if (s0 == 1 && s == 1)
            return 59;
          if (s0 == 1 && kp == 1 && s == 2)
            return 60;
          return 61;
        }
        return 99;
      },
      {"in_fiducial", "neutrino_pdg", "interaction_ccnc", "mc_n_strange",
       "count_kaon_plus", "count_kaon_minus", "count_kaon_zero",
       "count_lambda", "count_sigma_plus", "count_sigma_zero",
       "count_sigma_minus"});

  auto exclusive_alias_frame =
      exclusive_channel_frame.Define("exclusive_strange_channels", "excl_channel");
  return exclusive_alias_frame;
}

ROOT::RDF::RNode TruthClassifier::assign_channel_definitions(
    ROOT::RDF::RNode df) const {
  auto channel_frame = df.Define(
      "channel_definitions",
      [](bool fv, int nu, int ccnc, int s, int npi, int np, int npi0,
         int ngamma) {
        if (!fv) {
          if (nu == 0)
            return 1;
          return 2;
        }
        if (ccnc == 1)
          return 14;
        if (ccnc == 0 && s > 0) {
          if (s == 1)
            return 15;
          return 16;
        }
        if (std::abs(nu) == 12 && ccnc == 0)
          return 17;
        if (std::abs(nu) == 14 && ccnc == 0) {
          if (npi == 0 && np > 0)
            return 10;
          if (npi == 1 && npi0 == 0)
            return 11;
          if (npi0 > 0 || ngamma >= 2)
            return 12;
          if (npi > 1)
            return 13;
          return 18;
        }
        return 99;
      },
      {"in_fiducial", "neutrino_pdg", "interaction_ccnc", "mc_n_strange",
       "mc_n_pion", "mc_n_proton", "count_pi_zero", "count_gamma"});

  auto signal_frame = channel_frame.Define(
      "is_truth_signal", [](int ch) { return ch == 15 || ch == 16; },
      {"channel_definitions"});

  auto pure_signal_frame = signal_frame.Define(
      "pure_slice_signal",
      [](bool is_sig, float purity, float completeness) {
        return is_sig && purity > 0.5f && completeness > 0.1f;
      },
      {"is_truth_signal", "neutrino_purity_from_pfp",
       "neutrino_completeness_from_pfp"});

  return pure_signal_frame;
}

} // namespace faint
