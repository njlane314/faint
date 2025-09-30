#include "faint/SampleSet.h"

#include <cassert>

#include <nlohmann/json.hpp>

#include "faint/Log.h"
#include "faint/MuonSelector.h"

namespace faint {

SampleSet::SampleSet(const RunReader& runs, VariableRegistry variables,
                     const std::string& beam, std::vector<std::string> periods,
                     const std::string& ntuple_dir, bool blind)
    : runs_(runs),
      variables_(std::move(variables)),
      ntuple_dir_(ntuple_dir),
      beam_(beam),
      periods_(std::move(periods)),
      blind_(blind),
      total_pot_(0.0),
      total_triggers_(0) {
  build();
}

SampleSet::Map& SampleSet::frames() noexcept { return samples_; }

double SampleSet::total_pot() const noexcept { return total_pot_; }

long SampleSet::total_triggers() const noexcept { return total_triggers_; }

const std::string& SampleSet::beam() const noexcept { return beam_; }

const std::vector<std::string>& SampleSet::periods() const noexcept {
  return periods_;
}

const Run* SampleSet::run_for(const SampleKey& sk) const {
  auto it = run_cache_.find(sk);
  return it != run_cache_.end() ? it->second : nullptr;
}

void SampleSet::snapshot(const std::string& filter, const std::string& out,
                         const std::vector<std::string>& cols) const {
  snapshot_impl(filter, out, cols);
}

void SampleSet::snapshot(const Selection& selection, const std::string& out,
                         const std::vector<std::string>& cols) const {
  snapshot_impl(selection.str(), out, cols);
}

void SampleSet::print_branches() {
  log::debug("SampleSet::print_branches",
             "Available branches in loaded samples:");
  for (auto& [sample_key, sample_def] : samples_) {
    log::debug("SampleSet::print_branches", "--- Sample:", sample_key.str(),
               "---");
    auto branches = sample_def.nominal().GetColumnNames();
    for (const auto& branch : branches) {
      log::debug("SampleSet::print_branches", "  - ", branch);
    }
  }
}

namespace {

bool has_external_sample(const Run& run) {
  for (const auto& sample : run.samples) {
    if (sample.value("sample_type", "") == "ext") return true;
  }
  return false;
}

}  // namespace

void SampleSet::build() {
  std::vector<const Run*> to_process;
  const auto& all_runs = runs_.all();

#ifndef NDEBUG
  long expected_total_triggers = 0;
#endif

  for (auto& p : periods_) {
    const auto& rc = runs_.get(beam_, p);
    total_pot_ += rc.nominal_pot;
    total_triggers_ += rc.nominal_triggers;
    to_process.push_back(&rc);

#ifndef NDEBUG
    expected_total_triggers += rc.nominal_triggers;
#endif

    for (const auto& [key, candidate] : all_runs) {
      if (&candidate == &rc) continue;
      if (candidate.run_period != p) continue;
      if (!has_external_sample(candidate)) continue;

      const auto& er = candidate;
      total_pot_ += er.nominal_pot;
      total_triggers_ += er.nominal_triggers;
      to_process.push_back(&er);

#ifndef NDEBUG
      expected_total_triggers += er.nominal_triggers;
#endif
    }
  }

#ifndef NDEBUG
  assert(total_triggers_ == expected_total_triggers);
#endif

  for (const Run* rc : to_process) add_run(*rc);
}

void SampleSet::add_run(const Run& rc) {
  processors_.reserve(processors_.size() + rc.samples.size());
  for (auto& s : rc.samples) {
    if (s.contains("active") && !s.at("active").get<bool>()) {
      log::info("SampleSet::add_run", "Skipping inactive sample: ",
                s.at("sample_key").get<std::string>());
      continue;
    }

    auto pipeline = build_pipeline(s);
    processors_.push_back(std::move(pipeline));

    auto& proc = *processors_.back();
    Sample sample{s, rc.samples, ntuple_dir_, variables_, proc};
    SampleKey key = sample.key();

    run_cache_.emplace(key, &rc);
    samples_.emplace(std::move(key), std::move(sample));
  }
}

  std::unique_ptr<EventProcessor> SampleSet::build_pipeline(
    const nlohmann::json& sample) {
  auto weighter =
      std::make_unique<Weighter>(sample, total_pot_, total_triggers_);
  auto preselection = std::make_unique<PreSelection>();
  auto muon_selector = std::make_unique<MuonSelector>();
  auto truth_classifier = std::make_unique<TruthClassifier>();

  muon_selector->chain_processor(std::move(truth_classifier));
  preselection->chain_processor(std::move(muon_selector));
  weighter->chain_processor(std::move(preselection));

  return weighter;
}

void SampleSet::snapshot_impl(const std::string& filter, const std::string& out,
                              const std::vector<std::string>& cols) const {
  bool first = true;
  ROOT::RDF::RSnapshotOptions opts;
  for (auto const& [key, sample] : samples_) {
    auto df = sample.nominal();
    if (!filter.empty()) df = df.Filter(filter);
    opts.fMode = first ? "RECREATE" : "UPDATE";
    df.Snapshot(key.c_str(), out, cols, opts);
    first = false;
  }
}

}  // namespace faint
