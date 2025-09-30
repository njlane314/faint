#include "faint/Samples.h"

#include <cassert>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <utility>

#include <nlohmann/json.hpp>

#include "faint/Selection.h"
#include "faint/proc/MuonSelector.h"

namespace faint {
namespace {

ROOT::RDF::RNode open_frame(const std::string& base_dir, const std::string& rel,
                            EventProcessor& processor, sample::Origin origin) {
  auto path = base_dir + "/" + rel;
  ROOT::RDataFrame df("nuselection/EventSelectionFilter", path);
  return processor.process(df, origin);
}

ROOT::RDF::RNode filter_truth(ROOT::RDF::RNode df, const std::string& truth) {
  return truth.empty() ? df : df.Filter(truth);
}

ROOT::RDF::RNode exclude_truth(ROOT::RDF::RNode df,
                               const std::vector<std::string>& keys,
                               const nlohmann::json& all) {
  for (const auto& k : keys) {
    bool found = false;
    for (const auto& s : all) {
      if (s.at("sample_key").get<std::string>() == k) {
        auto filter_str = s.contains("truth")
                              ? s.at("truth").get<std::string>()
                              : s.value("truth_filter", std::string{});
        if (!filter_str.empty()) {
          df = df.Filter("!(" + filter_str + ")");
          found = true;
          break;
        }
      }
    }
    if (!found) {
      std::clog << "Sample::exclude_truth: Exclusion k not found or missing truth: "
                << k << '\n';
    }
  }
  return df;
}

bool has_external_sample(const Run& run) {
  for (const auto& sample : run.samples) {
    if (sample.value("sample_type", "") == "ext") return true;
  }
  return false;
}

}  // namespace

Sample::Sample(const nlohmann::json& j, const nlohmann::json& all,
               const std::string& base_dir, const VariableRegistry& vars,
               EventProcessor& processor)
    : key_{j.at("sample_key").get<std::string>()},
      origin_{[&]() {
        auto ts = j.at("sample_type").get<std::string>();
        if (ts == "mc") return sample::Origin::kMonteCarlo;
        if (ts == "data") return sample::Origin::kData;
        if (ts == "ext") return sample::Origin::kExternal;
        if (ts == "dirt") return sample::Origin::kDirt;
        return sample::Origin::kUnknown;
      }()},
      path_{j.value("relative_path", "")},
      truth_{[&]() {
        if (j.contains("truth")) return j.at("truth").get<std::string>();
        if (j.contains("truth_filter"))
          return j.at("truth_filter").get<std::string>();
        return std::string{};
      }()},
      exclude_{j.value("exclusion_truth_filters", std::vector<std::string>{})},
      pot_{j.value("pot", 0.0)},
      triggers_{j.value("triggers", 0L)},
      nominal_node_{build(base_dir, vars, processor, path_, all)} {
  if (j.contains("detector_variations")) {
    for (auto& dv : j.at("detector_variations")) {
      sample::Variation dvt =
          parse_variation(dv.at("variation_type").get<std::string>());
      variation_paths_[dvt] = dv.at("relative_path").get<std::string>();
    }
  }
  validate(base_dir);
  if (origin_ == sample::Origin::kMonteCarlo) {
    for (auto& [variation, rel_path] : variation_paths_) {
      variations_.emplace(variation,
                          build(base_dir, vars, processor, rel_path, all));
    }
  }
}

void Sample::validate(const std::string& base_dir) const {
  if (key_.str().empty()) {
    throw std::runtime_error("Sample::validate: empty key_");
  }
  if (origin_ == sample::Origin::kUnknown) {
    throw std::runtime_error("Sample::validate: unknown origin for " +
                             key_.str());
  }
  if ((origin_ == sample::Origin::kMonteCarlo || origin_ == sample::Origin::kDirt) &&
      pot_ <= 0) {
    std::ostringstream msg;
    msg << "Sample::validate: invalid pot_ for MC/Dirt " << key_.str()
        << " (pot=" << pot_ << ")";
    throw std::runtime_error(msg.str());
  }
  if (origin_ == sample::Origin::kData && triggers_ <= 0) {
    std::ostringstream msg;
    msg << "Sample::validate: invalid triggers_ for Data " << key_.str()
        << " (triggers=" << triggers_ << ")";
    throw std::runtime_error(msg.str());
  }
  if (origin_ != sample::Origin::kData && path_.empty()) {
    throw std::runtime_error("Sample::validate: missing path for " +
                             key_.str());
  }

  if (!path_.empty()) {
    auto p = std::filesystem::path(base_dir) / path_;
    if (!std::filesystem::exists(p)) {
      throw std::runtime_error("Sample::validate: missing file " + p.string());
    }
  }

  for (auto& [variation, rel_path] : variation_paths_) {
    auto vp = std::filesystem::path(base_dir) / rel_path;
    if (!std::filesystem::exists(vp)) {
      throw std::runtime_error("Sample::validate: missing variation " +
                               rel_path);
    }
  }
}

sample::Variation Sample::parse_variation(const std::string& s) const {
  if (s == "cv") return sample::Variation::kCV;
  if (s == "lyatt") return sample::Variation::kLYAttenuation;
  if (s == "lydown") return sample::Variation::kLYDown;
  if (s == "lyray") return sample::Variation::kLYRayleigh;
  if (s == "recomb2") return sample::Variation::kRecomb2;
  if (s == "sce") return sample::Variation::kSCE;
  if (s == "wiremodx") return sample::Variation::kWireModX;
  if (s == "wiremodyz") return sample::Variation::kWireModYZ;
  if (s == "wiremodanglexz") return sample::Variation::kWireModAngleXZ;
  if (s == "wiremodangleyz") return sample::Variation::kWireModAngleYZ;
  throw std::runtime_error("Sample::parse_variation: invalid detvar_type: " + s);
}

ROOT::RDF::RNode Sample::build(const std::string& base_dir,
                               [[maybe_unused]] const VariableRegistry& vars,
                               EventProcessor& processor,
                               const std::string& rel,
                               const nlohmann::json& all) {
  auto df = open_frame(base_dir, rel, processor, origin_);
  df = filter_truth(df, truth_);
  df = exclude_truth(df, exclude_, all);
  return df;
}

SampleSet::SampleSet(const RunReader& runs, VariableRegistry variables,
                     const std::string& beam, std::vector<std::string> periods,
                     const std::string& ntuple_dir)
    : runs_(runs),
      variables_(std::move(variables)),
      ntuple_dir_(ntuple_dir),
      beam_(beam),
      periods_(std::move(periods)),
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

const Run* SampleSet::run_for(const sample::Key& sk) const {
  auto it = run_cache_.find(sk);
  return it != run_cache_.end() ? it->second : nullptr;
}

void SampleSet::snapshot(const std::string& filter, const std::string& out,
                         const std::vector<std::string>& cols) const {
  snapshot_impl(filter, out, cols);
}

void SampleSet::snapshot(const selection::Selection& selection, const std::string& out,
                         const std::vector<std::string>& cols) const {
  snapshot_impl(selection.str(), out, cols);
}

void SampleSet::print_branches() {
#ifndef NDEBUG
  std::clog << "SampleSet::print_branches: Available branches in loaded samples:\n";
  for (auto& [sample_key, sample_def] : samples_) {
    std::clog << "SampleSet::print_branches: --- Sample: " << sample_key.str()
              << " ---\n";
    auto branches = sample_def.nominal().GetColumnNames();
    for (const auto& branch : branches) {
      std::clog << "SampleSet::print_branches:   - " << branch << '\n';
    }
  }
#endif
}

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
      std::clog << "SampleSet::add_run: Skipping inactive sample: "
                << s.at("sample_key").get<std::string>() << '\n';
      continue;
    }

    auto pipeline = build_pipeline(s);
    processors_.push_back(std::move(pipeline));

    auto& proc = *processors_.back();
    Sample sample{s, rc.samples, ntuple_dir_, variables_, proc};
    sample::Key key = sample.key();

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
