#include "faint/Dataset.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

namespace faint {
namespace dataset {

std::string run_config_path() {
    const char* env = gSystem->Getenv("FAINT_RUN_CONFIG");
    if (env) return env;
    std::string cwd = gSystem->pwd();
    return cwd + "/data/samples.json";
}

std::string ntuple_directory() {
    return ntuple_directory(run_config_path());
}

std::string ntuple_directory(const std::string& run_config_json) {
    const auto& config_path = run_config_json;
    std::ifstream input(config_path);
    if (!input.is_open()) {
        throw std::runtime_error("Could not open run configuration: " + config_path);
    }

    try {
        auto data = nlohmann::json::parse(input);
        if (data.contains("samples")) {
            data = data.at("samples");
        }

        if (!data.contains("ntupledir")) {
            throw std::runtime_error("Run configuration missing 'ntupledir' entry");
        }

        auto dir = data.at("ntupledir").get<std::string>();
        if (dir.empty()) {
            throw std::runtime_error("Run configuration has empty 'ntupledir'");
        }

        return dir;
    } catch (const nlohmann::json::exception& e) {
        throw std::runtime_error(std::string{"Failed to parse run configuration: "} + e.what());
    }
}

Dataset Dataset::open(const std::string& run_config_json, Options opt, Variables vars) {
    Dataset dataset;
    dataset.runs_ = RunReader::from_file(run_config_json);
    dataset.vars_ = std::move(vars);
    dataset.opt_ = std::move(opt);
    dataset.set_ = std::make_unique<SampleSet>(
        dataset.runs_, dataset.vars_, dataset.opt_.beam, dataset.opt_.periods,
        dataset.opt_.ntuple_dir, dataset.opt_.blind);
    dataset.build_dataset_cache();
    return dataset;
}

std::vector<std::string> Dataset::sample_keys(
    SampleOrigin origin_filter) const {
    std::vector<std::string> out;
    out.reserve(datasets_.size());
    for (auto const& [key, variations] : datasets_) {
        if (origin_filter == SampleOrigin::kUnknown ||
            variations.nominal.origin == origin_filter) {
            out.push_back(key.str());
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

ROOT::RDF::RNode Dataset::df(std::string_view sample_key,
                              SampleVariation v) const {
    const Variations* variations = find_dataset(sample_key);
    if (!variations) {
        throw std::runtime_error(std::string("Sample not found: ") + std::string(sample_key));
    }
    if (v == SampleVariation::kCV) return variations->nominal.dataframe();
    auto it = variations->variations.find(v);
    return (it != variations->variations.end()) ? it->second.dataframe()
                                                : variations->nominal.dataframe();
}

ROOT::RDF::RNode Dataset::final(std::string_view key,
                                 SampleVariation v) const {
    return df(key, v).Filter(sel::Final);
}

ROOT::RDF::RNode Dataset::quality(std::string_view key,
                                   SampleVariation v) const {
    return df(key, v).Filter(sel::Quality);
}

void Dataset::snapshot_where(const std::string& filter,
                              const std::string& out_file,
                              const std::vector<std::string>& columns) const {
    set().snapshot(filter, out_file, columns);
}

void Dataset::snapshot_final(const std::string& out_file,
                              const std::vector<std::string>& columns) const {
    snapshot_where(sel::Final, out_file, columns);
}

double Dataset::pot() const noexcept { return set().total_pot(); }

long Dataset::triggers() const noexcept { return set().total_triggers(); }

const std::string& Dataset::beam() const noexcept { return opt_.beam; }

const std::vector<std::string>& Dataset::periods() const noexcept {
    return opt_.periods;
}

const SampleSet& Dataset::samples() const { return set(); }

const RunReader& Dataset::runs() const { return runs_; }

const SampleSet& Dataset::set() const { return *set_; }

void Dataset::build_dataset_cache() {
    datasets_.clear();
    auto& frames = const_cast<SampleSet&>(set()).frames();
    for (auto& [key, sample] : frames) {
        Variations variations{
            make_entry(sample, SampleVariation::kCV, sample.nominal()), {}
        };
        for (auto const& [variation, node] : sample.variations()) {
            variations.variations.emplace(variation, make_entry(sample, variation, node));
        }
        datasets_.emplace(key, std::move(variations));
    }
}

const Dataset::Variations* Dataset::find_dataset(std::string_view key) const {
    const SampleKey sample_key{std::string(key)};
    auto it = datasets_.find(sample_key);
    return it != datasets_.end() ? &it->second : nullptr;
}

Dataset::Entry Dataset::make_entry(const Sample& sample, SampleVariation variation,
                                   ROOT::RDF::RNode node) {
    SampleRole role = SampleRole::kNominal;
    if (variation != SampleVariation::kCV) {
        role = SampleRole::kVariation;
    } else if (sample.origin() == SampleOrigin::kData) {
        role = SampleRole::kData;
    }

    return Entry{sample.origin(), role, std::move(node)};
}

} // namespace dataset
} // namespace faint
