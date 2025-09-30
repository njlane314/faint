#include "faint/proc/PreSelection.h"

#include "faint/Selections.h"

#include "ROOT/RVec.hxx"

namespace faint {

using sample::SampleOrigin;

ROOT::RDF::RNode PreSelection::process(ROOT::RDF::RNode df,
                                       SampleOrigin origin) const {
  ROOT::RDF::RNode node = df;

  if (!node.HasColumn("num_slices") && node.HasColumn("nslice"))
    node = node.Alias("num_slices", "nslice");

  if (!node.HasColumn("optical_filter_pe_beam")) {
    if (node.HasColumn("_opfilter_pe_beam"))
      node = node.Alias("optical_filter_pe_beam", "_opfilter_pe_beam");
    else if (node.HasColumn("opfilter_pe_beam"))
      node = node.Alias("optical_filter_pe_beam", "opfilter_pe_beam");
    else
      node = node.Define("optical_filter_pe_beam", []() { return 0.f; });
  }

  if (!node.HasColumn("optical_filter_pe_veto")) {
    if (node.HasColumn("_opfilter_pe_veto"))
      node = node.Alias("optical_filter_pe_veto", "_opfilter_pe_veto");
    else if (node.HasColumn("opfilter_pe_veto"))
      node = node.Alias("optical_filter_pe_veto", "opfilter_pe_veto");
    else
      node = node.Define("optical_filter_pe_veto", []() { return 0.f; });
  }

  if (!node.HasColumn("reco_neutrino_vertex_sce_x") &&
      node.HasColumn("reco_nu_vtx_sce_x"))
    node = node.Alias("reco_neutrino_vertex_sce_x", "reco_nu_vtx_sce_x");
  if (!node.HasColumn("reco_neutrino_vertex_sce_y") &&
      node.HasColumn("reco_nu_vtx_sce_y"))
    node = node.Alias("reco_neutrino_vertex_sce_y", "reco_nu_vtx_sce_y");
  if (!node.HasColumn("reco_neutrino_vertex_sce_z") &&
      node.HasColumn("reco_nu_vtx_sce_z"))
    node = node.Alias("reco_neutrino_vertex_sce_z", "reco_nu_vtx_sce_z");

  node = node.Define(
      "in_reco_fiducial",
      [](float x, float y, float z) {
        return selection::in_reco_fiducial_volume(x, y, z);
      },
      {"reco_neutrino_vertex_sce_x", "reco_neutrino_vertex_sce_y",
       "reco_neutrino_vertex_sce_z"});

  if (!node.HasColumn("n_pfps_gen2")) {
    node = node.Define(
        "n_pfps_gen2",
        [](const ROOT::RVec<unsigned> &gens) {
          return ROOT::VecOps::Sum(gens == 2u);
        },
        {"pfp_generations"});
  }
  if (!node.HasColumn("n_pfps_gen3")) {
    node = node.Define(
        "n_pfps_gen3",
        [](const ROOT::RVec<unsigned> &gens) {
          return ROOT::VecOps::Sum(gens == 3u);
        },
        {"pfp_generations"});
  }

  if (origin == SampleOrigin::kMonteCarlo) {
    if (node.HasColumn("software_trigger_pre_ext")) {
      node = node.Define(
          "software_trigger",
          [](unsigned run, int pre, int post) {
            return run < 16880 ? pre > 0 : post > 0;
          },
          {"run", "software_trigger_pre_ext", "software_trigger_post_ext"});
    } else if (node.HasColumn("software_trigger_pre")) {
      node = node.Define(
          "software_trigger",
          [](unsigned run, int pre, int post) {
            return run < 16880 ? pre > 0 : post > 0;
          },
          {"run", "software_trigger_pre", "software_trigger_post"});
    } else if (node.HasColumn("software_trigger")) {
      node = node.Redefine("software_trigger", "software_trigger != 0");
    } else {
      node = node.Define("software_trigger", []() { return true; });
    }
  } else {
    if (node.HasColumn("software_trigger"))
      node = node.Redefine("software_trigger", "software_trigger != 0");
    else
      node = node.Define("software_trigger", []() { return true; });
  }

  node = node.Define(
      selection::column::kPassPre,
      [origin](float pe_beam, float pe_veto, bool swtrig) {
        return selection::passes_pre_selection(origin, pe_beam, pe_veto,
                                               swtrig);
      },
      {"optical_filter_pe_beam", "optical_filter_pe_veto",
       "software_trigger"});

  node = node.Define(
      selection::column::kPassFlash,
      [](int nslices, float topo, int n_gen2) {
        return selection::passes_flash_selection(nslices, topo, n_gen2);
      },
      {"num_slices", "topological_score", "n_pfps_gen2"});

  node = node.Define(
      selection::column::kPassFiducial,
      [](float x, float y, float z) {
        return selection::in_reco_fiducial_volume(x, y, z);
      },
      {"reco_neutrino_vertex_sce_x", "reco_neutrino_vertex_sce_y",
       "reco_neutrino_vertex_sce_z"});

  if (!node.HasColumn("n_muons_tot"))
    node = node.Define("n_muons_tot", []() { return 0UL; });
  node = node.Define(
      selection::column::kPassMuon,
      [](unsigned long n_muons) {
        return selection::passes_muon_selection(n_muons);
      },
      {"n_muons_tot"});

  node = node.Define(
      selection::column::kPassTopology,
      [](float contained, float cluster) {
        return selection::passes_topology_selection(contained, cluster);
      },
      {"contained_fraction", "slice_cluster_fraction"});

  node = node.Define(
      selection::column::kPassFinal,
      [](bool pre, bool flash, bool fv, bool mu, bool topo) {
        return selection::passes_final_selection(pre, flash, fv, mu, topo);
      },
      {selection::column::kPassPre, selection::column::kPassFlash,
       selection::column::kPassFiducial, selection::column::kPassMuon,
       selection::column::kPassTopology});

  node = node.Define(
      selection::column::kQualityEvent,
      [](bool pre, bool flash, bool fv, bool topo) {
        return selection::is_quality_event(pre, flash, fv, topo);
      },
      {selection::column::kPassPre, selection::column::kPassFlash,
       selection::column::kPassFiducial, selection::column::kPassTopology});

  return next_ ? next_->process(node, origin) : node;
}

} // namespace faint
