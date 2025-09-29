#include "faint/SampleSet.h"

#include <nlohmann/json.hpp>

#include "faint/MuonSelector.h"
#include "faint/utils/Logger.h"

namespace faint {

SampleSet::SampleSet(const RunCatalog& runs, VariableRegistry variables,
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

void SampleSet::snapshot(const SelectionQuery& query, const std::string& out,
                         const std::vector<std::string>& cols) const {
  snapshot_impl(query.str(), out, cols);
}

void SampleSet::print_branches() {
  log::debug("SampleSet::print_branches",
             "Available branches in loaded samples:");
  for (auto& [sample_key, sample_def] : samples_) {
    log::debug("SampleSet::print_branches", "--- Sample:", sample_key.str(),
               "---");
    auto branches = sample_def.nominal_node_.GetColumnNames();
    for (const auto& branch : branches) {
      log::debug("SampleSet::print_branches", "  - ", branch);
    }
  }
}

void SampleSet::build() {
  const std::string ext_beam{"numi_ext"};
  std::vector<const Run*> to_process;

  for (auto& p : periods_) {
    const auto& rc = runs_.get(beam_, p);
    total_pot_ += rc.nominal_pot;
    total_triggers_ += rc.nominal_triggers;
    to_process.push_back(&rc);

    auto key = ext_beam + ":" + p;
    if (runs_.all().count(key)) {
      const auto& er = runs_.get(ext_beam, p);
      total_pot_ += er.nominal_pot;
      total_triggers_ += er.nominal_triggers;
      to_process.push_back(&er);
    }
  }

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
    Samples sample{s, rc.samples, ntuple_dir_, variables_, proc};

    run_cache_.emplace(sample.sample_key_, &rc);
    samples_.emplace(sample.sample_key_, std::move(sample));
  }
}

std::unique_ptr<IEventProcessor> SampleSet::build_pipeline(
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
    auto df = sample.nominal_node_;
    if (!filter.empty()) df = df.Filter(filter);
    opts.fMode = first ? "RECREATE" : "UPDATE";
    df.Snapshot(key.c_str(), out, cols, opts);
    first = false;
  }
}

}  // namespace faint
