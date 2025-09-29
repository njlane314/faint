#include "faint/Campaign.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

namespace faint {
namespace campaign {

std::string run_config_path() {
    const char* env = gSystem->Getenv("FAINT_RUN_CONFIG");
    if (env) return env;
    std::string cwd = gSystem->pwd();
    return cwd + "/data/samples.json";
}

std::string ntuple_directory() {
    const char* env = gSystem->Getenv("FAINT_NTUPLES");
    if (!env) throw std::runtime_error("Set FAINT_NTUPLES to the directory containing faint ntuples");
    return env;
}

Campaign Campaign::open(const std::string& run_config_json, Options opt, Variables vars) {
    Campaign s;
    s.runs_ = RunCatalog::from_file(run_config_json);
    s.vars_ = std::move(vars);
    s.opt_ = std::move(opt);
    s.set_ = std::make_unique<SampleSet>(
        s.runs_, s.vars_, s.opt_.beam, s.opt_.periods, s.opt_.ntuple_dir, s.opt_.blind);
    return s;
}

std::vector<std::string> Campaign::sample_keys(
    SampleOrigin origin_filter) const {
  std::vector<std::string> out;
  const auto& frames = const_cast<SampleSet&>(set()).frames();
  for (auto const& kv : frames) {
    const auto& key = kv.first;
    const auto& sample = kv.second;
    if (origin_filter == SampleOrigin::kUnknown ||
        sample.origin() == origin_filter) {
      out.push_back(key.str());
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

ROOT::RDF::RNode Campaign::df(std::string_view sample_key,
                              SampleVariation v) const {
  const Sample* sample = find_sample(sample_key);
  if (!sample) {
    throw std::runtime_error(std::string("Sample not found: ") + std::string(sample_key));
  }
  if (v == SampleVariation::kCV) return sample->nominal();
  const auto& variations = sample->variations();
  auto it = variations.find(v);
  return (it != variations.end()) ? it->second : sample->nominal();
}

ROOT::RDF::RNode Campaign::final(std::string_view key,
                                 SampleVariation v) const {
  return df(key, v).Filter(sel::Final);
}

ROOT::RDF::RNode Campaign::quality(std::string_view key,
                                   SampleVariation v) const {
  return df(key, v).Filter(sel::Quality);
}

void Campaign::snapshot_where(const std::string& filter,
                              const std::string& out_file,
                              const std::vector<std::string>& columns) const {
    set().snapshot(filter, out_file, columns);
}

void Campaign::snapshot_final(const std::string& out_file,
                              const std::vector<std::string>& columns) const {
    snapshot_where(sel::Final, out_file, columns);
}

double Campaign::pot() const noexcept { return set().total_pot(); }

long Campaign::triggers() const noexcept { return set().total_triggers(); }

const std::string& Campaign::beam() const noexcept { return opt_.beam; }

const std::vector<std::string>& Campaign::periods() const noexcept {
    return opt_.periods;
}

const SampleSet& Campaign::samples() const { return set(); }

const RunCatalog& Campaign::runs() const { return runs_; }

const SampleSet& Campaign::set() const { return *set_; }

const Sample* Campaign::find_sample(std::string_view key) const {
    const auto& frames = const_cast<SampleSet&>(set()).frames();
    const SampleKey sample_key{std::string(key)};
    auto it = frames.find(sample_key);
    return it != frames.end() ? &it->second : nullptr;
}

} // namespace campaign
} // namespace faint
