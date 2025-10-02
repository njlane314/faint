#pragma once
#include <ROOT/RDataFrame.hxx>
#include <RtypesCore.h>
#include <iostream>
#include <string>
#include <vector>

#include "rarexsec/Volume.hh"
#include "rarexsec/Hub.hh"

namespace rarexsec {
namespace selection {

namespace cuts {
inline constexpr float min_beam_pe = 0.f;
inline constexpr float max_veto_pe = 20.f;
inline constexpr int required_slices = 1;
inline constexpr float min_topological_score = 0.06f;
inline constexpr int min_generation2_pfps = 2;
inline constexpr float min_contained_fraction = 0.7f;
inline constexpr float min_cluster_fraction = 0.5f;
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
                                  ? (pe_beam > cuts::min_beam_pe &&
                                     pe_veto < cuts::max_veto_pe)
                                  : true;
    return dataset_gate && software_trigger;
}

inline bool passes_flash_selection(int num_slices, float topological_score,
                                   int generation2_pfps) {
    return num_slices == cuts::required_slices &&
           topological_score > cuts::min_topological_score &&
           generation2_pfps >= cuts::min_generation2_pfps;
}

inline bool in_reco_fiducial_volume(float x, float y, float z) {
    return fiducial::is_in_reco_volume(x, y, z);
}

inline bool passes_muon_selection(std::size_t n_muons) {
    return n_muons > 0;
}

inline bool passes_topology_selection(float contained_fraction,
                                      float cluster_fraction) {
    return contained_fraction >= cuts::min_contained_fraction &&
           cluster_fraction >= cuts::min_cluster_fraction;
}

inline bool passes_muon_track_selection(float s, float llr, float l, float d,
                                        unsigned g) {
    return s > cuts::min_score && llr > cuts::min_llr && l > cuts::min_length &&
           d < cuts::max_distance && g == cuts::required_generation;
}

inline bool passes_final_selection(bool pre_ok, bool flash_ok, bool fiducial_ok,
                                   bool muon_ok, bool topology_ok) {
    return pre_ok && flash_ok && fiducial_ok && muon_ok && topology_ok;
}

inline bool is_quality_event(bool pre_ok, bool flash_ok, bool fiducial_ok, bool topology_ok) {
    return pre_ok && flash_ok && fiducial_ok && topology_ok;
}

enum class Preset {
    All,
    FiducialOnly,
    MuonOnly,
    Baseline,
    PreOnly,
    FlashOnly,
    TopologyOnly,
    Final
};

inline ROOT::RDF::RNode apply(ROOT::RDF::RNode node, Preset p, const rarexsec::Entry& rec) {
    switch (p) {
        case Preset::All:
            return node;
        case Preset::FiducialOnly:
            return node.Filter([](bool fv){ return fv; },
                               {"in_reco_fiducial"});
        case Preset::MuonOnly:
            return node.Filter([](bool mu){ return mu; },
                               {"has_muon"});
        case Preset::Baseline:
            return node.Filter([](bool fv, bool mu){ return fv && mu; },
                               {"in_reco_fiducial", "has_muon"});
        case Preset::PreOnly:
            return node.Filter([k = rec.kind](float pe_beam, float pe_veto, bool sw){
                                   return passes_pre_selection(k, pe_beam, pe_veto, sw);
                               },
                               {"pe_beam", "pe_veto", "software_trigger"});
        case Preset::FlashOnly:
            return node.Filter([](int ns, float topo, int n2g){
                                   return passes_flash_selection(ns, topo, n2g);
                               },
                               {"num_slices", "topological_score", "generation2_pfps"});
        case Preset::TopologyOnly:
            return node.Filter([](float cf, float cl){
                                   return passes_topology_selection(cf, cl);
                               },
                               {"contained_fraction", "cluster_fraction"});
        case Preset::Final:
        default:
            return node.Filter([k = rec.kind](float pe_beam, float pe_veto, bool sw,
                                              int ns, float topo, int n2g,
                                              bool fv, bool mu,
                                              float cf, float cl){
                                   const bool pre_ok = passes_pre_selection(k, pe_beam, pe_veto, sw);
                                   const bool flash_ok = passes_flash_selection(ns, topo, n2g);
                                   const bool topo_ok = passes_topology_selection(cf, cl);
                                   return passes_final_selection(pre_ok, flash_ok, fv, mu, topo_ok);
                               },
                               {"pe_beam", "pe_veto", "software_trigger",
                                "num_slices", "topological_score", "generation2_pfps",
                                "in_reco_fiducial", "has_muon",
                                "contained_fraction", "cluster_fraction"});
    }
}

}
}
