#ifndef SAMPLE_DEFINITION_H
#define SAMPLE_DEFINITION_H

#include <filesystem>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "nlohmann/json.hpp"

#include "faint/EventProcessor.h"
#include "faint/Logger.h"
#include "faint/Types.h"
#include "faint/Variables.h"

namespace faint {

inline ROOT::RDF::RNode open_frame(const std::string& base_dir,
                                   const std::string& rel,
                                   EventProcessor& processor,
                                   SampleOrigin origin) {
  auto path = base_dir + "/" + rel;
  ROOT::RDataFrame df("nuselection/EventSelectionFilter", path);
  return processor.process(df, origin);
}

inline ROOT::RDF::RNode filter_truth(ROOT::RDF::RNode df, const std::string& truth) {
  return truth.empty() ? std::move(df) : df.Filter(truth);
}

inline ROOT::RDF::RNode exclude_truth(ROOT::RDF::RNode df,
                                      const std::vector<std::string>& keys,
                                      const nlohmann::json& all) {
  ROOT::RDF::RNode node = std::move(df);
  for (const auto& key : keys) {
    bool found = false;
    for (const auto& sample_json : all) {
      if (sample_json.at("sample_key").get<std::string>() == key) {
        if (sample_json.contains("truth")) {
          auto filter_str = sample_json.at("truth").get<std::string>();
          node = node.Filter("!(" + filter_str + ")");
        }
        found = true;
        break;
      }
    }
    if (!found) {
      log::warn("Sample::exclude_truth", "Exclusion key not found or missing truth:", key);
    }
  }
  return node;
}

class Sample {
 public:
  SampleKey key_{};
  SampleOrigin origin_{SampleOrigin::kUnknown};

  std::string path_;
  std::string truth_;
  std::vector<std::string> exclude_;

  double pot_{0.0};
  long triggers_{0};

  ROOT::RDF::RNode node_;
  std::map<SampleVariation, ROOT::RDF::RNode> variations_;

  Sample(const nlohmann::json& sample_json,
         const nlohmann::json& all_samples,
         const std::string& base_dir,
         const VariableRegistry& vars,
         EventProcessor& processor)
      : key_{sample_json.at("sample_key").get<std::string>()},
        origin_{[&]() {
          auto type = sample_json.at("sample_type").get<std::string>();
          if (type == "mc") return SampleOrigin::kMonteCarlo;
          if (type == "data") return SampleOrigin::kData;
          if (type == "ext") return SampleOrigin::kExternal;
          if (type == "dirt") return SampleOrigin::kDirt;
          return SampleOrigin::kUnknown;
        }()},
        path_{sample_json.value("relative_path", "")},
        truth_{sample_json.value("truth", sample_json.value("truth_filter", ""))},
        exclude_{sample_json.value("exclusion_truth_filters", std::vector<std::string>{})},
        pot_{sample_json.value("pot", 0.0)},
        triggers_{sample_json.value("triggers", 0L)},
        node_{build(base_dir, vars, processor, path_, all_samples)} {
    if (sample_json.contains("detector_variations")) {
      for (const auto& detvar : sample_json.at("detector_variations")) {
        auto variation = parse_variation(detvar.at("variation_type").get<std::string>());
        variation_paths_[variation] = detvar.at("relative_path").get<std::string>();
      }
    }

    validate(base_dir);

    if (origin_ == SampleOrigin::kMonteCarlo) {
      for (auto& [variation, rel_path] : variation_paths_) {
        variations_.emplace(variation, build(base_dir, vars, processor, rel_path, all_samples));
      }
    }
  }

  void validate(const std::string& base_dir) const {
    if (key_.empty()) log::fatal("Sample::validate", "Empty sample key");
    if (origin_ == SampleOrigin::kUnknown) log::fatal("Sample::validate", "Unknown sample origin for", key_.str());

    if ((origin_ == SampleOrigin::kMonteCarlo || origin_ == SampleOrigin::kDirt) && pot_ <= 0.0) {
      log::fatal("Sample::validate", "Non-positive POT for", key_.str());
    }

    if (origin_ == SampleOrigin::kData && triggers_ <= 0L) {
      log::fatal("Sample::validate", "Non-positive trigger count for", key_.str());
    }

    if (origin_ != SampleOrigin::kData && path_.empty()) {
      log::fatal("Sample::validate", "Missing relative path for", key_.str());
    }

    if (!path_.empty()) {
      auto full_path = std::filesystem::path(base_dir) / path_;
      if (!std::filesystem::exists(full_path)) {
        log::fatal("Sample::validate", "Missing sample file", full_path.string());
      }
    }

    for (const auto& [variation, rel_path] : variation_paths_) {
      auto variation_full_path = std::filesystem::path(base_dir) / rel_path;
      if (!std::filesystem::exists(variation_full_path)) {
        log::fatal("Sample::validate", "Missing variation file", rel_path);
      }
    }
  }

 private:
  std::map<SampleVariation, std::string> variation_paths_;

  SampleVariation parse_variation(const std::string& value) const {
    if (value == "cv") return SampleVariation::kCV;
    if (value == "lyatt") return SampleVariation::kLYAttenuation;
    if (value == "lydown") return SampleVariation::kLYDown;
    if (value == "lyray") return SampleVariation::kLYRayleigh;
    if (value == "recomb2") return SampleVariation::kRecomb2;
    if (value == "sce") return SampleVariation::kSCE;
    if (value == "wiremodx") return SampleVariation::kWireModX;
    if (value == "wiremodyz") return SampleVariation::kWireModYZ;
    if (value == "wiremodanglexz") return SampleVariation::kWireModAngleXZ;
    if (value == "wiremodangleyz") return SampleVariation::kWireModAngleYZ;
    log::fatal("Sample::parse_variation", "Unsupported variation type", value);
    return SampleVariation::kUnknown;
  }

  ROOT::RDF::RNode build(const std::string& base_dir,
                         const VariableRegistry& vars,
                         EventProcessor& processor,
                         const std::string& rel,
                         const nlohmann::json& all) {
    auto node = open_frame(base_dir, rel, processor, origin_);
    node = filter_truth(std::move(node), truth_);
    node = exclude_truth(std::move(node), exclude_, all);
    (void)vars;
    return node;
  }
};

}  // namespace faint

#endif
