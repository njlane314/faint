#pragma once

#include <cstddef>
#include <string>

#include "rarexsec/Volume.hh"
#include "rarexsec/Hub.hh"

namespace rarexsec {
namespace selection {

namespace column {
inline constexpr const char* pass_pre = "pass_pre";
inline constexpr const char* pass_flash = "pass_flash";
inline constexpr const char* pass_fiducial = "pass_fv";
inline constexpr const char* pass_muon = "pass_mu";
inline constexpr const char* pass_topology = "pass_topo";
inline constexpr const char* pass_final = "pass_final";
inline constexpr const char* quality_event = "quality_event";
inline constexpr const char* nominal_weight = "nominal_event_weight";
}

namespace pre {
inline constexpr float min_beam_pe = 0.f;
inline constexpr float max_veto_pe = 20.f;
}

namespace flash {
inline constexpr int required_slices = 1;
inline constexpr float min_topological_score = 0.06f;
inline constexpr int min_generation2_pfps = 2;
}

namespace topology {
inline constexpr float min_contained_fraction = 0.7f;
inline constexpr float min_cluster_fraction = 0.5f;
}

namespace muon_track {
inline constexpr float min_score = 0.5f;
inline constexpr float min_llr = 0.2f;
inline constexpr float min_length = 10.0f;
inline constexpr float max_distance = 4.0f;
inline constexpr unsigned required_generation = 2u;
}

inline bool passes_pre_selection(sample::origin origin, float pe_beam,
                                 float pe_veto, bool software_trigger) {
    const bool requires_dataset_gate =
        (origin == sample::origin::beam || origin == sample::origin::strangeness ||
         origin == sample::origin::dirt);
    const bool dataset_gate = requires_dataset_gate
                                  ? (pe_beam > pre::min_beam_pe &&
                                     pe_veto < pre::max_veto_pe)
                                  : true;
    return dataset_gate && software_trigger;
}

inline bool passes_flash_selection(int num_slices, float topological_score,
                                   int generation2_pfps) {
    return num_slices == flash::required_slices &&
           topological_score > flash::min_topological_score &&
           generation2_pfps >= flash::min_generation2_pfps;
}

inline bool in_reco_fiducial_volume(float x, float y, float z) {
    return fiducial::is_in_reco_volume(x, y, z);
}

inline bool passes_muon_selection(std::size_t n_muons) {
    return n_muons > 0;
}

inline bool passes_topology_selection(float contained_fraction,
                                      float cluster_fraction) {
    return contained_fraction >= topology::min_contained_fraction &&
           cluster_fraction >= topology::min_cluster_fraction;
}

inline bool passes_muon_track_selection(float s, float llr, float l, float d, unsigned g) {
    return s > muon_track::min_score && llr > muon_track::min_llr && l > muon_track::min_length && d < muon_track::max_distance && g == muon_track::required_generation;
}

inline bool passes_final_selection(bool pre, bool flash, bool fiducial,
                                   bool muon, bool topology) {
    return pre && flash && fiducial && muon && topology;
}

inline bool is_quality_event(bool pre, bool flash, bool fiducial, bool topology) {
    return pre && flash && fiducial && topology;
}

class Selection {
public:
    Selection();
    explicit Selection(std::string expression);

    const std::string& str() const noexcept;
    bool empty() const noexcept;

private:
    std::string expression_;
};

}
}
