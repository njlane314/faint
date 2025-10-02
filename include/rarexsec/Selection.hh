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

inline const char* to_string(Preset p) {
    switch (p) {
        case Preset::All: return "All";
        case Preset::FiducialOnly: return "FiducialOnly";
        case Preset::MuonOnly: return "MuonOnly";
        case Preset::Baseline: return "Baseline";
        case Preset::PreOnly: return "PreOnly";
        case Preset::FlashOnly: return "FlashOnly";
        case Preset::TopologyOnly: return "TopologyOnly";
        case Preset::Final: return "Final";
        default: return "Unknown";
    }
}

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

struct DetVarResult {
    std::string tag;
    ULong64_t entries = 0;
};

struct SampleResult {
    const rarexsec::Entry* entry = nullptr;
    ULong64_t entries = 0;
    double weighted = 0.0;
    std::vector<DetVarResult> detvars;
};

struct Summary {
    std::vector<SampleResult> samples;
    ULong64_t total_entries = 0;
    double total_weighted = 0.0;
    double total_pot_nom = 0.0;
    double total_pot_eqv = 0.0;
    double total_trig_nom = 0.0;
    double total_trig_eqv = 0.0;
};

inline Summary evaluate(const std::vector<const rarexsec::Entry*>& samples,
                        Preset preset,
                        const char* weight_col = "w_nominal") {
    Summary out;
    out.samples.reserve(samples.size());
    for (const auto* e : samples) {
        if (!e) continue;
        auto r_sel = apply(e->rnode(), preset, *e);
        auto n_ptr = r_sel.Count();
        auto w_ptr = r_sel.Sum<double>(weight_col);
        const auto n_val = n_ptr.GetValue();
        const auto w_val = w_ptr.GetValue();
        SampleResult row;
        row.entry = e;
        row.entries = n_val;
        row.weighted = w_val;
        for (const auto& kv : e->detvars) {
            const auto& tag = kv.first;
            const auto& dat = kv.second;
            if (!dat.node) continue;
            auto r_dv = apply(dat.rnode(), preset, *e);
            const auto dv_n = r_dv.Count().GetValue();
            row.detvars.push_back({tag, dv_n});
        }
        out.samples.push_back(std::move(row));
        out.total_entries += n_val;
        out.total_weighted += w_val;
        out.total_pot_nom += e->pot_nom;
        out.total_pot_eqv += (e->pot_eqv > 0.0 ? e->pot_eqv : e->pot_nom);
        out.total_trig_nom += e->trig_nom;
        out.total_trig_eqv += (e->trig_eqv > 0.0 ? e->trig_eqv : e->trig_nom);
    }
    return out;
}

inline Summary evaluate(const rarexsec::Hub& hub,
                        const std::string& beamline,
                        const std::vector<std::string>& periods,
                        Preset preset,
                        const char* weight_col = "w_nominal") {
    auto sims = hub.simulation(beamline, periods);
    return evaluate(sims, preset, weight_col);
}

inline const char* origin_to_string(rarexsec::sample::origin k) {
    switch (k) {
        case rarexsec::sample::origin::data: return "data";
        case rarexsec::sample::origin::beam: return "beam";
        case rarexsec::sample::origin::strangeness: return "strangeness";
        case rarexsec::sample::origin::ext: return "ext";
        case rarexsec::sample::origin::dirt: return "dirt";
        default: return "unknown";
    }
}

inline void print(const Summary& S, std::ostream& os = std::cout) {
    for (const auto& row : S.samples) {
        os << "Sample kind '" << origin_to_string(row.entry->kind)
           << "' from file " << row.entry->file << "\n";
        os << "  Final selection entries: " << row.entries
           << " | weighted: " << row.weighted << "\n";
        for (const auto& dv : row.detvars) {
            os << "  Detector variation '" << dv.tag
               << "' entries: " << dv.entries << "\n";
        }
    }
    os << "Total POT (nominal): " << S.total_pot_nom << "\n"
       << "Total POT (equivalent): " << S.total_pot_eqv << "\n"
       << "Total triggers (nominal): " << S.total_trig_nom << "\n"
       << "Total triggers (equivalent): " << S.total_trig_eqv << "\n";
}

} 
}