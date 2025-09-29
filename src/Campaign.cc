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
    s.runs_ = RunCatalog::fromFile(run_config_json);
    s.vars_ = std::move(vars);
    s.opt_ = std::move(opt);
    s.set_ = std::make_unique<SampleSet>(
        s.runs_, s.vars_, s.opt_.beam, s.opt_.periods, s.opt_.ntuple_dir, s.opt_.blind);
    return s;
}

std::vector<std::string> Campaign::sample_keys(Origin origin_filter) const {
    std::vector<std::string> out;
    const auto& frames = const_cast<SampleSet&>(set()).frames();
    for (auto const& kv : frames) {
        const auto& key = kv.first;
        const auto& sample = kv.second;
        if (origin_filter == Origin::kUnknown || sample.origin_ == origin_filter) {
            out.push_back(key.str());
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

ROOT::RDF::RNode Campaign::df(std::string_view sample_key, Variation v) const {
    const Sample* sample = find_sample(sample_key);
    if (!sample) {
        throw std::runtime_error(std::string("Sample not found: ") + std::string(sample_key));
    }
    if (v == Variation::kCV) return sample->node_;
    auto it = sample->variations_.find(v);
    return (it != sample->variations_.end()) ? it->second : sample->node_;
}

ROOT::RDF::RNode Campaign::final(std::string_view key, Variation v) const {
    return df(key, v).Filter(sel::Final);
}

ROOT::RDF::RNode Campaign::quality(std::string_view key, Variation v) const {
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

const Sample* Campaign::find_sample(std::string_view key) const {
    const auto& frames = const_cast<SampleSet&>(set()).frames();
    for (auto const& kv : frames) {
        if (kv.first.str() == key) return &kv.second;
    }
    return nullptr;
}

} // namespace campaign
} // namespace faint
