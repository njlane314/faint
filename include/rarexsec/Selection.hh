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
                                   const bool requires_dataset_gate =
                                       (k == sample::origin::beam ||
                                        k == sample::origin::strangeness ||
                                        k == sample::origin::dirt);
                                   const bool dataset_gate = requires_dataset_gate
                                                                ? (pe_beam > cuts::min_beam_pe &&
                                                                   pe_veto < cuts::max_veto_pe)
                                                                : true;
                                   return dataset_gate && sw;
                               },
                               {"pe_beam", "pe_veto", "software_trigger"});
        case Preset::FlashOnly:
            return node.Filter([](int ns, float topo, int n2g){
                                   return ns == cuts::required_slices &&
                                          topo > cuts::min_topological_score &&
                                          n2g >= cuts::min_generation2_pfps;
                               },
                               {"num_slices", "topological_score", "generation2_pfps"});
        case Preset::TopologyOnly:
            return node.Filter([](float cf, float cl){
                                   return cf >= cuts::min_contained_fraction &&
                                          cl >= cuts::min_cluster_fraction;
                               },
                               {"contained_fraction", "cluster_fraction"});
        case Preset::Final:
        default:
            return node.Filter([k = rec.kind](float pe_beam, float pe_veto, bool sw,
                                              int ns, float topo, int n2g,
                                              bool fv, bool mu,
                                              float cf, float cl){
                                   const bool requires_dataset_gate =
                                       (k == sample::origin::beam ||
                                        k == sample::origin::strangeness ||
                                        k == sample::origin::dirt);
                                   const bool dataset_gate = requires_dataset_gate
                                                                ? (pe_beam > cuts::min_beam_pe &&
                                                                   pe_veto < cuts::max_veto_pe)
                                                                : true;
                                   const bool pre_ok = dataset_gate && sw;
                                   const bool flash_ok = (ns == cuts::required_slices &&
                                                          topo > cuts::min_topological_score &&
                                                          n2g >= cuts::min_generation2_pfps);
                                   const bool topo_ok = (cf >= cuts::min_contained_fraction &&
                                                         cl >= cuts::min_cluster_fraction);
                                   return pre_ok && flash_ok && fv && mu && topo_ok;
                               },
                               {"pe_beam", "pe_veto", "software_trigger",
                                "num_slices", "topological_score", "generation2_pfps",
                                "in_reco_fiducial", "has_muon",
                                "contained_fraction", "cluster_fraction"});
    }
}

}
}
