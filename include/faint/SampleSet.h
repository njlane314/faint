#ifndef ANALYSIS_SAMPLE_SET_H
#define ANALYSIS_SAMPLE_SET_H

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ROOT/RDataFrame.hxx"

#include "faint/EventProcessor.h"
#include "faint/Logger.h"
#include "faint/MuonSelector.h"
#include "faint/PreSelection.h"
#include "faint/Run.h"
#include "faint/RunCatalog.h"
#include "faint/Sample.h"
#include "faint/TruthClassifier.h"
#include "faint/Variables.h"
#include "faint/Weighter.h"

namespace faint {

class SampleSet {
 public:
  using Map = std::map<SampleKey, Sample>;

  SampleSet(const RunCatalog& runs,
            VariableRegistry variables,
            std::string beam,
            std::vector<std::string> periods,
            std::string ntuple_dir,
            bool blind = true)
      : runs_(runs),
        variables_(std::move(variables)),
        ntuple_dir_(std::move(ntuple_dir)),
        beam_(std::move(beam)),
        periods_(std::move(periods)),
        blind_(blind) {
    build();
  }

  Map& frames() noexcept { return samples_; }
  const Map& frames() const noexcept { return samples_; }

  double total_pot() const noexcept { return total_pot_; }
  long total_triggers() const noexcept { return total_triggers_; }

  const std::string& beam() const noexcept { return beam_; }
  const std::vector<std::string>& periods() const noexcept { return periods_; }

  const Run* run_for(const SampleKey& key) const {
    auto it = run_cache_.find(key);
    return it != run_cache_.end() ? it->second : nullptr;
  }

  void snapshot(const std::string& filter,
                const std::string& out,
                const std::vector<std::string>& columns = {}) const {
    bool first = true;
    ROOT::RDF::RSnapshotOptions opts;
    for (const auto& [key, sample] : samples_) {
      auto node = sample.node_;
      if (!filter.empty()) {
        node = node.Filter(filter);
      }
      opts.fMode = first ? "RECREATE" : "UPDATE";
      node.Snapshot(key.str().c_str(), out, columns, opts);
      first = false;
    }
  }

  void snapshot_final(const std::string& out,
                      const std::vector<std::string>& columns = {}) const {
    snapshot("pass_final", out, columns);
  }

  void print_branches() const {
    log::debug("SampleSet::print_branches", "Available branches in loaded samples:");
    for (const auto& [sample_key, sample] : samples_) {
      log::debug("SampleSet::print_branches", "--- Sample:", sample_key.str(), "---");
      auto branches = sample.node_.GetColumnNames();
      for (const auto& branch : branches) {
        log::debug("SampleSet::print_branches", "  - ", branch);
      }
    }
  }

 private:
  const RunCatalog& runs_;
  VariableRegistry variables_;
  std::string ntuple_dir_;
  std::string beam_;
  std::vector<std::string> periods_;
  bool blind_;

  double total_pot_{0.0};
  long total_triggers_{0};

  Map samples_;
  std::vector<std::unique_ptr<EventProcessor>> processors_;
  std::unordered_map<SampleKey, const Run*> run_cache_;

  void build() {
    for (const auto& period : periods_) {
      const Run& run = runs_.get(beam_, period);
      add_run(run);

      auto ext_beam = infer_external_beam(beam_);
      if (!ext_beam.empty()) {
        auto key = ext_beam + ":" + period;
        if (runs_.all().count(key)) {
          add_run(runs_.get(ext_beam, period));
        }
      }
    }
  }

  static std::string infer_external_beam(const std::string& beam) {
    auto dash = beam.find('-');
    if (dash == std::string::npos) {
      return beam + "-ext";
    }
    return beam.substr(0, dash) + "-ext";
  }

  void add_run(const Run& run) {
    total_pot_ += run.nominal_pot;
    total_triggers_ += run.nominal_triggers;

    processors_.reserve(processors_.size() + run.samples.size());
    for (const auto& sample_json : run.samples) {
      if (sample_json.contains("active") && !sample_json.at("active").get<bool>()) {
        log::info("SampleSet::add_run", "Skipping inactive sample:", sample_json.at("sample_key").get<std::string>());
        continue;
      }

      auto pipeline = chain(std::make_unique<Weighter>(sample_json, total_pot_, total_triggers_),
                            std::make_unique<PreSelection>(),
                            std::make_unique<MuonSelector>(),
                            std::make_unique<TruthClassifier>());

      processors_.push_back(std::move(pipeline));
      auto& processor = *processors_.back();

      Sample sample{sample_json, run.samples, ntuple_dir_, variables_, processor};
      run_cache_.emplace(sample.key_, &run);
      samples_.emplace(sample.key_, std::move(sample));
    }
  }

  template <typename Head, typename... Tail>
  std::unique_ptr<EventProcessor> chain(std::unique_ptr<Head> head, std::unique_ptr<Tail>... tail) {
    if constexpr (sizeof...(tail) == 0) {
      return head;
    } else {
      auto next = chain(std::move(tail)...);
      head->chain_processor(std::move(next));
      return head;
    }
  }
};

}  // namespace faint

#endif
