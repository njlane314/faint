#include "faint/Sample.h"

#include <filesystem>
#include <utility>

#include "faint/Log.h"

namespace faint {
namespace {

ROOT::RDF::RNode open_frame(const std::string& base_dir, const std::string& rel,
                            EventProcessor& processor, SampleOrigin origin) {
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
        if (s.contains("truth")) {
          auto filter_str = s.at("truth").get<std::string>();
          df = df.Filter("!(" + filter_str + ")");
          found = true;
          break;
        }
      }
    }
    if (!found)
      log::warn("Sample::exclude_truth", "Exclusion k not found or missing truth:", k);
  }
  return df;
}

}  // namespace

Sample::Sample(const nlohmann::json& j, const nlohmann::json& all,
               const std::string& base_dir, const VariableRegistry& vars,
               EventProcessor& processor)
    : key_{j.at("sample_key").get<std::string>()},
      origin_{[&]() {
        auto ts = j.at("sample_type").get<std::string>();
        if (ts == "mc") return SampleOrigin::kMonteCarlo;
        if (ts == "data") return SampleOrigin::kData;
        if (ts == "ext") return SampleOrigin::kExternal;
        if (ts == "dirt") return SampleOrigin::kDirt;
        return SampleOrigin::kUnknown;
      }()},
      path_{j.value("relative_path", "")},
      truth_{j.value("truth", "")},
      exclude_{j.value("exclusion_truth_filters", std::vector<std::string>{})},
      pot_{j.value("pot", 0.0)},
      triggers_{j.value("triggers", 0L)},
      nominal_node_{build(base_dir, vars, processor, path_, all)} {
  if (j.contains("detector_variations")) {
    for (auto& dv : j.at("detector_variations")) {
      SampleVariation dvt =
          parse_variation(dv.at("variation_type").get<std::string>());
      variation_paths_[dvt] = dv.at("relative_path").get<std::string>();
    }
  }
  validate(base_dir);
  if (origin_ == SampleOrigin::kMonteCarlo) {
    for (auto& [variation, rel_path] : variation_paths_) {
      variations_.emplace(variation,
                          build(base_dir, vars, processor, rel_path, all));
    }
  }
}

void Sample::validate(const std::string& base_dir) const {
  if (key_.str().empty())
    log::fatal("Sample::validate", "empty key_");
  if (origin_ == SampleOrigin::kUnknown)
    log::fatal("Sample::validate", "unknown origin for", key_.str());
  if ((origin_ == SampleOrigin::kMonteCarlo || origin_ == SampleOrigin::kDirt) &&
      pot_ <= 0)
    log::fatal("Sample::validate", "invalid pot_ for MC/Dirt", key_.str());
  if (origin_ == SampleOrigin::kData && triggers_ <= 0)
    log::fatal("Sample::validate", "invalid triggers_ for Data", key_.str());
  if (origin_ != SampleOrigin::kData && path_.empty())
    log::fatal("Sample::validate", "missing path for", key_.str());

  if (!path_.empty()) {
    auto p = std::filesystem::path(base_dir) / path_;
    if (!std::filesystem::exists(p))
      log::fatal("Sample::validate", "missing file", p.string());
  }

  for (auto& [variation, rel_path] : variation_paths_) {
    auto vp = std::filesystem::path(base_dir) / rel_path;
    if (!std::filesystem::exists(vp))
      log::fatal("Sample::validate", "missing variation", rel_path);
  }
}

SampleVariation Sample::parse_variation(const std::string& s) const {
  if (s == "cv") return SampleVariation::kCV;
  if (s == "lyatt") return SampleVariation::kLYAttenuation;
  if (s == "lydown") return SampleVariation::kLYDown;
  if (s == "lyray") return SampleVariation::kLYRayleigh;
  if (s == "recomb2") return SampleVariation::kRecomb2;
  if (s == "sce") return SampleVariation::kSCE;
  if (s == "wiremodx") return SampleVariation::kWireModX;
  if (s == "wiremodyz") return SampleVariation::kWireModYZ;
  if (s == "wiremodanglexz") return SampleVariation::kWireModAngleXZ;
  if (s == "wiremodangleyz") return SampleVariation::kWireModAngleYZ;
  log::fatal("Sample::parse_variation", "invalid detvar_type:", s);
  return SampleVariation::kUnknown;
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

}  // namespace faint
