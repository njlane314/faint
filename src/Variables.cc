#include "faint/Variables.h"

#include "faint/Samples.h"

#include <unordered_set>
#include <utility>

namespace faint {

using sample::SampleOrigin;

const Variables::KnobVariations& Variables::knob_var() {
  static const KnobVariations m = {
      {"RPA", {"knobRPAup", "knobRPAdn"}},
      {"CCMEC", {"knobCCMECup", "knobCCMECdn"}},
      {"AxFFCCQE", {"knobAxFFCCQEup", "knobAxFFCCQEdn"}},
      {"VecFFCCQE", {"knobVecFFCCQEup", "knobVecFFCCQEdn"}},
      {"DecayAngMEC", {"knobDecayAngMECup", "knobDecayAngMECdn"}},
      {"ThetaDelta2Npi", {"knobThetaDelta2Npiup", "knobThetaDelta2Npidn"}},
      {"ThetaDelta2NRad", {"knobThetaDelta2NRadup", "knobThetaDelta2NRaddn"}},
      {"NormCCCOH", {"knobNormCCCOHup", "knobNormCCCOHdn"}},
      {"NormNCCOH", {"knobNormNCCOHup", "knobNormNCCOHdn"}},
      {"xsr_scc_Fv3", {"knobxsr_scc_Fv3up", "knobxsr_scc_Fv3dn"}},
      {"xsr_scc_Fa3", {"knobxsr_scc_Fa3up", "knobxsr_scc_Fa3dn"}}};
  return m;
}

const Variables::MultiUniverseVars& Variables::multi_uni_var() {
  static const MultiUniverseVars m = {{"weightsGenie", 500},
                                      {"weightsFlux", 500},
                                      {"weightsReint", 500},
                                      {"weightsPPFX", 500}};
  return m;
}

const std::string& Variables::single_knob_var() {
  static const std::string s = "RootinoFix";
  return s;
}

std::vector<std::string> Variables::event_var(SampleOrigin origin) {
  auto vars = collect_base_vars();
  if (origin == SampleOrigin::kMonteCarlo || origin == SampleOrigin::kDirt)
    add_mc_vars(vars);
  return std::vector<std::string>(vars.begin(), vars.end());
}

std::unordered_set<std::string> Variables::collect_base_vars() {
  std::unordered_set<std::string> vars(base_var().begin(), base_var().end());
  vars.insert(reco_var().begin(), reco_var().end());
  vars.insert(track_var().begin(), track_var().end());
  vars.insert(proc_evt_var().begin(), proc_evt_var().end());
  vars.insert(image_var().begin(), image_var().end());
  vars.insert(flash_var().begin(), flash_var().end());
  vars.insert(energy_var().begin(), energy_var().end());
  vars.insert(slice_var().begin(), slice_var().end());
  return vars;
}

void Variables::add_mc_vars(std::unordered_set<std::string>& vars) {
  vars.insert(truth_var().begin(), truth_var().end());
  for (const auto& kv : knob_var()) {
    vars.insert(kv.second.first);
    vars.insert(kv.second.second);
  }
  for (const auto& kv : multi_uni_var())
    vars.insert(kv.first);
  vars.insert(single_knob_var());
}

const std::vector<std::string>& Variables::base_var() {
  static const std::vector<std::string> v = {"run", "sub", "evt"};
  return v;
}

const std::vector<std::string>& Variables::truth_var() {
  static const std::vector<std::string> v = {
      "neutrino_pdg",        "interaction_ccnc",      "interaction_mode",
      "interaction_type",    "neutrino_energy",       "neutrino_theta",
      "neutrino_pt",         "target_nucleus_pdg",    "hit_nucleon_pdg",
      "kinematic_W",         "kinematic_X",           "kinematic_Y",
      "kinematic_Q_squared", "neutrino_momentum_x",   "neutrino_momentum_y",
      "neutrino_momentum_z", "neutrino_vertex_x",     "neutrino_vertex_y",
      "neutrino_vertex_z",   "neutrino_vertex_wire_u","neutrino_vertex_wire_v",
      "neutrino_vertex_wire_w","neutrino_vertex_time", "neutrino_sce_vertex_x",
      "neutrino_sce_vertex_y","neutrino_sce_vertex_z", "lepton_energy",
      "true_neutrino_momentum_x","true_neutrino_momentum_y",
      "true_neutrino_momentum_z", "flux_path_length",  "flux_parent_pdg",
      "flux_hadron_pdg",     "flux_decay_mode",       "flux_decay_vtx_x",
      "flux_decay_vtx_y",    "flux_decay_vtx_z",      "flux_decay_mom_x",
      "flux_decay_mom_y",    "flux_decay_mom_z",      "numi_baseline",
      "numi_off_axis_angle", "bnb_baseline",          "bnb_off_axis_angle",
      "is_vertex_in_fiducial","count_mu_minus",       "count_mu_plus",
      "count_e_minus",       "count_e_plus",          "count_pi_zero",
      "count_pi_plus",       "count_pi_minus",        "count_kaon_plus",
      "count_kaon_minus",    "count_kaon_zero",       "count_proton",
      "count_neutron",       "count_gamma",           "count_lambda",
      "count_sigma_plus",    "count_sigma_zero",      "count_sigma_minus",
      "mc_particle_pdg",     "mc_particle_trackid",   "mc_particle_energy",
      "mc_elastic_scatters", "mc_inelastic_scatters", "mc_momentum_x",
      "mc_momentum_y",       "mc_momentum_z",         "mc_end_momentum",
      "mc_start_vertex_x",   "mc_start_vertex_y",     "mc_start_vertex_z",
      "mc_end_vertex_x",     "mc_end_vertex_y",       "mc_end_vertex_z",
      "mc_particle_final_state","mc_completeness",    "mc_purity",
      "mc_daughter_pdg",     "mc_daughter_energy",    "mc_daughter_process_flat",
      "mc_daughter_process_idx","mc_daughter_mom_x",  "mc_daughter_mom_y",
      "mc_daughter_mom_z",   "mc_daughter_vtx_x",     "mc_daughter_vtx_y",
      "mc_daughter_vtx_z",   "mc_allchain_primary_index","mc_allchain_trackid",
      "mc_allchain_pdg",     "mc_allchain_energy",    "mc_allchain_elastic_scatters",
      "mc_allchain_inelastic_scatters","mc_allchain_momentum_x",
      "mc_allchain_momentum_y","mc_allchain_momentum_z","mc_allchain_end_momentum",
      "mc_allchain_start_vertex_x","mc_allchain_start_vertex_y",
      "mc_allchain_start_vertex_z", "mc_allchain_end_vertex_x",
      "mc_allchain_end_vertex_y",    "mc_allchain_end_vertex_z",
      "mc_allchain_parent_trackid",  "mc_allchain_process",
      "mc_allchain_final_state",     "mc_allchain_completeness",
      "mc_allchain_purity",          "true_transverse_momentum",
      "true_visible_transverse_momentum", "true_total_momentum",
      "true_visible_total_momentum", "true_visible_energy",
      "neutrino_completeness_from_pfp", "neutrino_purity_from_pfp",
      "backtracked_pdg_codes",        "blip_pdg"};
  return v;
}

const std::vector<std::string>& Variables::reco_var() {
  static const std::vector<std::string> v = {
      "reco_neutrino_vertex_sce_x", "reco_neutrino_vertex_sce_y",
      "reco_neutrino_vertex_sce_z", "num_slices",
      "slice_num_hits", "selection_pass", "slice_id",
      "optical_filter_pe_beam", "optical_filter_pe_veto",
      "num_pfps", "num_tracks", "num_showers", "event_total_hits",
      "crt_veto", "crt_hit_pe", "pfp_slice_indices",
      "backtracked_pdg_codes", "backtracked_energies",
      "backtracked_track_ids", "backtracked_purities",
      "backtracked_completenesses", "backtracked_overlay_purities",
      "backtracked_momentum_x", "backtracked_momentum_y",
      "backtracked_momentum_z", "backtracked_start_x",
      "backtracked_start_y", "backtracked_start_z", "backtracked_start_time",
      "backtracked_start_wire_U", "backtracked_start_wire_V",
      "backtracked_start_wire_Y", "backtracked_sce_start_x",
      "backtracked_sce_start_y", "backtracked_sce_start_z",
      "backtracked_sce_start_wire_U", "backtracked_sce_start_wire_V",
      "backtracked_sce_start_wire_Y", "software_trigger",
      "software_trigger_pre", "software_trigger_post",
      "software_trigger_pre_ext", "software_trigger_post_ext",
      "slice_pdg", "pfp_generations", "pfp_track_daughters",
      "pfp_shower_daughters", "pfp_num_descendents", "pfp_vertex_x",
      "pfp_vertex_y", "pfp_vertex_z", "track_shower_scores",
      "pfp_pdg_codes", "pfp_num_hits", "pfp_num_plane_hits_U",
      "pfp_num_plane_hits_V", "pfp_num_plane_hits_Y", "pfp_num_subclusters_U",
      "pfp_num_subclusters_V", "pfp_num_subclusters_Y",
      "pfp_max_subhit_fraction_U", "pfp_max_subhit_fraction_V",
      "pfp_max_subhit_fraction_Y", "total_hits_U", "total_hits_V",
      "total_hits_Y", "slice_topological_scores", "topological_score",
      "slice_cluster_fraction", "contained_fraction"};
  return v;
}

const std::vector<std::string>& Variables::image_var() {
  static const std::vector<std::string> v = {
      "reco_neutrino_vertex_x",  "reco_neutrino_vertex_y",
      "reco_neutrino_vertex_z",  "detector_image_u",
      "detector_image_v",        "detector_image_w",
      "semantic_image_u",        "semantic_image_v",
      "semantic_image_w",        "event_detector_image_u",
      "event_detector_image_v",  "event_detector_image_w",
      "event_semantic_image_u",  "event_semantic_image_v",
      "event_semantic_image_w",  "event_adc_u",
      "event_adc_v",             "event_adc_w",
      "slice_semantic_counts_u", "slice_semantic_counts_v",
      "slice_semantic_counts_w", "event_semantic_counts_u",
      "event_semantic_counts_v", "event_semantic_counts_w",
      "is_vtx_in_image_u",       "is_vtx_in_image_v",
      "is_vtx_in_image_w",       "inference_score"};
  return v;
}

const std::vector<std::string>& Variables::flash_var() {
  static const std::vector<std::string> v = {
      "t0",            "flash_match_score", "flash_total_pe",
      "flash_time",    "flash_z_center",    "flash_z_width",
      "slice_charge",  "slice_z_center",    "charge_light_ratio",
      "flash_slice_z_dist", "flash_pe_per_charge"};
  return v;
}

const std::vector<std::string>& Variables::energy_var() {
  static const std::vector<std::string> v = {
      "neutrino_energy_0", "neutrino_energy_1", "neutrino_energy_2",
      "slice_calo_energy_0", "slice_calo_energy_1", "slice_calo_energy_2"};
  return v;
}

const std::vector<std::string>& Variables::slice_var() {
  static const std::vector<std::string> v = {
      "original_event_neutrino_hits", "event_neutrino_hits",
      "event_muon_hits",              "event_electron_hits",
      "event_proton_hits",            "event_charged_pion_hits",
      "event_neutral_pion_hits",      "event_neutron_hits",
      "event_gamma_hits",             "event_other_hits",
      "event_charged_kaon_hits",      "event_neutral_kaon_hits",
      "event_lambda_hits",            "event_charged_sigma_hits",
      "event_sigma_zero_hits",        "event_cosmic_hits",
      "slice_neutrino_hits",          "slice_muon_hits",
      "slice_electron_hits",          "slice_proton_hits",
      "slice_charged_pion_hits",      "slice_neutral_pion_hits",
      "slice_neutron_hits",           "slice_gamma_hits",
      "slice_other_hits",             "slice_charged_kaon_hits",
      "slice_neutral_kaon_hits",      "slice_lambda_hits",
      "slice_charged_sigma_hits",     "slice_sigma_zero_hits",
      "slice_cosmic_hits",            "pfp_neutrino_hits",
      "pfp_muon_hits",                "pfp_electron_hits",
      "pfp_proton_hits",              "pfp_charged_pion_hits",
      "pfp_neutral_pion_hits",        "pfp_neutron_hits",
      "pfp_gamma_hits",               "pfp_other_hits",
      "pfp_charged_kaon_hits",        "pfp_neutral_kaon_hits",
      "pfp_lambda_hits",              "pfp_charged_sigma_hits",
      "pfp_sigma_zero_hits",          "pfp_cosmic_hits",
      "neutrino_completeness_from_pfp", "neutrino_purity_from_pfp"};
  return v;
}

const std::vector<std::string>& Variables::track_var() {
  static const std::vector<std::string> v = {
      "track_shower_scores",     "trk_llr_pid_v",
      "track_length",            "track_distance_to_vertex",
      "track_start_x",           "track_start_y",
      "track_start_z",           "track_end_x",
      "track_end_y",             "track_end_z",
      "track_theta",             "track_phi",
      "track_calo_energy_u",     "track_calo_energy_v",
      "track_calo_energy_y"};
  return v;
}

const std::vector<std::string>& Variables::proc_evt_var() {
  static const std::vector<std::string> v = {
      "in_reco_fiducial",         "reco_nu_vtx_sce_x",
      "reco_nu_vtx_sce_y",        "reco_nu_vtx_sce_z",
      "n_pfps_gen2",              "n_pfps_gen3",
      "quality_event",            "n_muons_tot",
      "has_muon",                  "muon_trk_score_v",
      "muon_trk_llr_pid_v",       "muon_trk_start_x_v",
      "muon_trk_start_y_v",       "muon_trk_start_z_v",
      "muon_trk_end_x_v",         "muon_trk_end_y_v",
      "muon_trk_end_z_v",         "muon_trk_length_v",
      "muon_trk_distance_v",      "muon_pfp_generation_v",
      "muon_trk_range_muon_mom_v","muon_track_costheta",
      "base_event_weight",        "nominal_event_weight",
      "in_fiducial",              "mc_n_strange",
      "mc_n_pion",                "mc_n_proton",
      "genie_int_mode",           "incl_channel",
      "excl_channel"};
  return v;
}

}  // namespace faint
