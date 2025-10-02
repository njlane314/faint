#pragma once
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <RtypesCore.h>
#include <iostream>
#include <cstddef>
#include <string>
#include <vector>

#include "rarexsec/Volume.hh"
#include "rarexsec/Hub.hh"

namespace rarexsec {
namespace selection {

inline constexpr float trigger_min_beam_pe = 0.f;
inline constexpr float trigger_max_veto_pe = 20.f;

inline constexpr int slice_required_count = 1;
inline constexpr float slice_min_topology_score = 0.06f;

inline constexpr float topology_min_contained_fraction = 0.7f;
inline constexpr float topology_min_cluster_fraction = 0.5f;

inline constexpr float muon_min_track_score = 0.5f;
inline constexpr float muon_min_llr = 0.2f;
inline constexpr float muon_min_track_length = 10.0f;
inline constexpr float muon_max_track_distance = 4.0f;
inline constexpr unsigned muon_required_generation = 2u;

enum class Preset {
    Empty,
    Trigger,
    Slice,
    Fiducial,
    Topology,
    Muon,
    InclusiveMuCC
};

inline ROOT::RDF::RNode apply(ROOT::RDF::RNode node, Preset p, const rarexsec::Entry& rec) {
    switch (p) {
        case Preset::Empty:
            return node;
        case Preset::Trigger:
            return node.Filter([k = rec.kind](float pe_beam, float pe_veto, bool sw){
                                   const bool requires_dataset_gate =
                                       (k == sample::origin::beam ||
                                        k == sample::origin::strangeness ||
                                        k == sample::origin::dirt);
                                   const bool dataset_gate = requires_dataset_gate
                                                                ? (pe_beam > trigger_min_beam_pe &&
                                                                   pe_veto < trigger_max_veto_pe)
                                                                : true;
                                   return dataset_gate && sw;
                               },
                               {"pe_beam", "pe_veto", "software_trigger"});
        case Preset::Slice:
            return node.Filter([](int ns, float topo){
                                   return ns == required_slices &&
                                          topo > min_topological_score;
                               },
                               {"num_slices", "topological_score"});
        case Preset::Fiducial:
            return node.Filter([](bool fv){ return fv; },
                               {"in_reco_fiducial"});
        case Preset::Topology:
            return node.Filter([](float cf, float cl){
                                   return cf >= topology_min_contained_fraction &&
                                          cl >= topology_min_cluster_fraction;
                               },
                               {"contained_fraction", "cluster_fraction"});
        case Preset::Muon:
            return node.Filter(
                [](const ROOT::RVec<float>& scores,
                   const ROOT::RVec<float>& llrs,
                   const ROOT::RVec<float>& lengths,
                   const ROOT::RVec<float>& distances,
                   const ROOT::RVec<unsigned>& generations) {
                    const auto n = scores.size();
                    for (std::size_t i = 0; i < n; ++i) {
                        const bool passes = scores[i] > muon_min_track_score &&
                                            llrs[i] > muon_min_llr &&
                                            lengths[i] > muon_min_track_length &&
                                            distances[i] < muon_max_track_distance &&
                                            generations[i] == muon_required_generation;
                        if (passes) {
                            return true;
                        }
                    }
                    return false;
                },
                {"track_shower_scores",
                 "trk_llr_pid_v",
                 "track_length",
                 "track_distance_to_vertex",
                 "pfp_generations"});
        case Preset::InclusiveMuCC:
        default: {
            auto filtered = apply(node, Preset::Trigger, rec);
            filtered = apply(filtered, Preset::Slice, rec);
            filtered = apply(filtered, Preset::Fiducial, rec);
            filtered = apply(filtered, Preset::Topology, rec);
            return apply(filtered, Preset::Muon, rec);
        }
    }
}

}
}
