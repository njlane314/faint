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

enum class Preset {
    Empty,
    Trigger,
    Slice, 
    FiducialOnly,
    MuonOnly,
    FlashOnly,
    TopologyOnly,
    Final
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
                                                                ? (pe_beam > min_beam_pe &&
                                                                   pe_veto < max_veto_pe)
                                                                : true;
                                   return dataset_gate && sw;
                               },
                               {"pe_beam", "pe_veto", "software_trigger"});
        // define a slice preset here
        case Preset::FiducialOnly:
            return node.Filter([](bool fv){ return fv; },
                               {"in_reco_fiducial"});
        case Preset::MuonOnly:
            return node.Filter(
                [](const ROOT::RVec<float>& scores,
                   const ROOT::RVec<float>& llrs,
                   const ROOT::RVec<float>& lengths,
                   const ROOT::RVec<float>& distances,
                   const ROOT::RVec<unsigned>& generations) {
                    const auto n = scores.size();
                    for (std::size_t i = 0; i < n; ++i) {
                        const bool passes = scores[i] > min_score &&
                                            llrs[i] > min_llr &&
                                            lengths[i] > min_length &&
                                            distances[i] < max_distance &&
                                            generations[i] == required_generation;
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

       
        case Preset::FlashOnly:
            return node.Filter([](int ns, float topo, int n2g){
                                   return ns == required_slices &&
                                          topo > min_topological_score &&
                                          n2g >= min_generation2_pfps;
                               },
                               {"num_slices", "topological_score", "generation2_pfps"});
        case Preset::TopologyOnly:
            return node.Filter([](float cf, float cl){
                                   return cf >= min_contained_fraction &&
                                          cl >= min_cluster_fraction;
                               },
                               {"contained_fraction", "cluster_fraction"});
        case Preset::Final:
        default: {
            auto filtered = apply(node, Preset::Trigger, rec);
            filtered = apply(filtered, Preset::FlashOnly, rec);
            filtered = apply(filtered, Preset::TopologyOnly, rec);
            return apply(filtered, Preset::Baseline, rec);
        }
    }
}

}
}
