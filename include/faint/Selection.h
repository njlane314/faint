#ifndef FAINT_SELECTION_H
#define FAINT_SELECTION_H

#include <cstddef>
#include <string>

#include <faint/FiducialVolume.h>
#include <faint/Types.h>

namespace faint {

class Selection {
 public:
  Selection();
  explicit Selection(std::string expression);

  const std::string& str() const noexcept;
  bool empty() const noexcept;

 private:
  std::string expression_;
};

namespace selection {

namespace column {
inline constexpr const char* kPassPre = "pass_pre";
inline constexpr const char* kPassFlash = "pass_flash";
inline constexpr const char* kPassFiducial = "pass_fv";
inline constexpr const char* kPassMuon = "pass_mu";
inline constexpr const char* kPassTopology = "pass_topo";
inline constexpr const char* kPassFinal = "pass_final";
inline constexpr const char* kQualityEvent = "quality_event";
}  // namespace column

struct PreCut {
  static constexpr float kMinBeamPE = 0.f;
  static constexpr float kMaxVetoPE = 20.f;
};

struct FlashCut {
  static constexpr int kRequiredSlices = 1;
  static constexpr float kMinTopologicalScore = 0.06f;
  static constexpr int kMinGeneration2PFPs = 2;
};

struct TopologyCut {
  static constexpr float kMinContainedFraction = 0.7f;
  static constexpr float kMinClusterFraction = 0.5f;
};

struct MuonTrackCut {
  static constexpr float kMinScore = 0.5f;
  static constexpr float kMinLLR = 0.2f;
  static constexpr float kMinLength = 10.0f;
  static constexpr float kMaxDistance = 4.0f;
  static constexpr unsigned kRequiredGeneration = 2u;
};

inline bool passes_pre_selection(SampleOrigin origin, float pe_beam,
                                 float pe_veto, bool software_trigger) {
  const bool requires_dataset_gate =
      (origin == SampleOrigin::kMonteCarlo || origin == SampleOrigin::kDirt);
  const bool dataset_gate = requires_dataset_gate
                                ? (pe_beam > PreCut::kMinBeamPE &&
                                   pe_veto < PreCut::kMaxVetoPE)
                                : true;
  return dataset_gate && software_trigger;
}

inline bool passes_flash_selection(int num_slices, float topological_score,
                                   int generation2_pfps) {
  return num_slices == FlashCut::kRequiredSlices &&
         topological_score > FlashCut::kMinTopologicalScore &&
         generation2_pfps >= FlashCut::kMinGeneration2PFPs;
}

inline bool in_reco_fiducial_volume(float x, float y, float z) {
  return fiducial::is_in_reco_volume(x, y, z);
}

inline bool passes_muon_selection(std::size_t n_muons) {
  return n_muons > 0;
}

inline bool passes_topology_selection(float contained_fraction,
                                      float cluster_fraction) {
  return contained_fraction >= TopologyCut::kMinContainedFraction &&
         cluster_fraction >= TopologyCut::kMinClusterFraction;
}

inline bool passes_muon_track_selection(float score, float llr, float length,
                                        float distance, unsigned generation,
                                        bool fid_start, bool fid_end) {
  return score > MuonTrackCut::kMinScore &&
         llr > MuonTrackCut::kMinLLR &&
         length > MuonTrackCut::kMinLength &&
         distance < MuonTrackCut::kMaxDistance &&
         generation == MuonTrackCut::kRequiredGeneration && fid_start &&
         fid_end;
}

inline bool passes_final_selection(bool pre, bool flash, bool fiducial,
                                   bool muon, bool topology) {
  return pre && flash && fiducial && muon && topology;
}

inline bool is_quality_event(bool pre, bool flash, bool fiducial,
                             bool topology) {
  return pre && flash && fiducial && topology;
}

}  // namespace selection

}  // namespace faint

#endif  // FAINT_SELECTION_H
