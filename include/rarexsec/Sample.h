#ifndef SAMPLE_DEFINITION_H
#define SAMPLE_DEFINITION_H

#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "nlohmann/json.hpp"

#include <rarexsec/utils/Logger.h>
#include <rarexsec/data/IEventProcessor.h>
#include <rarexsec/data/SampleTypes.h>
#include <rarexsec/data/VariableRegistry.h>

namespace analysis {

inline ROOT::RDF::RNode open_frame(const std::string &base_dir, 
                                   const std::string &rel,
                                   IEventProcessor &processor, 
                                   SampleOrigin origin) {
    auto path = base_dir + "/" + rel;
    ROOT::RDataFrame df("nuselection/EventSelectionFilter", path);
    return processor.process(df, origin);
}

inline ROOT::RDF::RNode filter_truth(ROOT::RDF::RNode df, const std::string &truth) {
    return truth.empty() ? df : df.Filter(truth);
}

inline ROOT::RDF::RNode exclude_truth(ROOT::RDF::RNode df, 
                                      const std::vector<std::string> &keys,
                                      const nlohmann::json &all) {
    for (const auto &k : keys) {
        bool found = false;
        for (const auto &s : all) {
            if (s.at("sample_key").get<std::string>() == k) {
                if (s.contains("truth")) {
                    auto filter_str = s.at("truth").get<std::string>();
                    df = df.Filter("!(" + filter_str + ")");
                    found = true;
                    break;
                }
            }
        }
        if (!found) log::warn("Sample::exclude_truth", "Exclusion k not found or missing truth:", k);
    }
    return df;
}

class Sample {
  public:
    SampleKey key_;
    SampleOrigin origin_;

    std::string path_;
    std::string truth_;
    std::vector<std::string> exclude_;

    double pot_{0.0};
    long triggers_{0};

    ROOT::RDF::RNode node_;
    std::map<SampleVariation, ROOT::RDF::RNode> variations_;

    Sample(const nlohmann::json &j, 
           const nlohmann::json &all, 
           const std::string &base_dir,
           const VariableRegistry &vars, 
           IEventProcessor &processor)
        : key_{j.at("sample_key").get<std::string>()},
          origin_{[&]() {
              auto ts = j.at("sample_type").get<std::string>();
              return (ts == "mc"     ? SampleOrigin::kMonteCarlo
                      : ts == "data" ? SampleOrigin::kData
                      : ts == "ext"  ? SampleOrigin::kExternal
                      : ts == "dirt" ? SampleOrigin::kDirt
                                      : SampleOrigin::kUnknown);
          }()},
          path_{j.value("relative_path", "")},
          truth_{j.value("truth", "")},
          exclude_{j.value("exclusion_truth_filters", std::vector<std::string>{})},
          pot_{j.value("pot", 0.0)},
          triggers_{j.value("triggers", 0L)},
          node_{build(base_dir, vars, processor, path_, all)} {
        if (j.contains("detector_variations")) {
            for (auto &dv : j.at("detector_variations")) {
                SampleVariation dvt = this->parse_variation(dv.at("variation_type").get<std::string>());
                variation_paths_[dvt] = dv.at("relative_path").get<std::string>();
            }
        }
        this->validate(base_dir);
        if (origin_ == SampleOrigin::kMonteCarlo) {
            for (auto &[v, path] : variation_paths_) {
                variations_.emplace(v, this->build(base_dir, vars, processor, path, all));
            }
        }
    }

    void validate(const std::string &base_dir) const {
        if (key_.str().empty()) log::fatal("Sample::validate", "empty key_");
        if (origin_ == SampleOrigin::kUnknown) log::fatal("Sample::validate", "unknown  for", key_.str());
        if ((origin_ == SampleOrigin::kMonteCarlo || origin_ == SampleOrigin::kDirt) && pot_ <= 0) log::fatal("Sample::validate", "invalid pot_ for MC/Dirt", key_.str());
        if (origin_ == SampleOrigin::kData && triggers_ <= 0) log::fatal("Sample::validate", "invalid triggers_ for Data", key_.str());
        if (origin_ != SampleOrigin::kData && path_.empty()) log::fatal("Sample::validate", "missing path for", key_.str());
        
        if (!path_.empty()) {
            auto p = std::filesystem::path(base_dir) / path_;
            if (!std::filesystem::exists(p))
                log::fatal("Sample::validate", "missing file", p.string());
        }

        for (auto &[v, rp] : variation_paths_) {
            auto vp = std::filesystem::path(base_dir) / rp;
            if (!std::filesystem::exists(vp)) log::fatal("Sample::validate", "missing variation", rp);
        }
    }

  private:
    std::map<SampleVariation, std::string> variation_paths_;

    SampleVariation parse_variation(const std::string &s) const {
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

    ROOT::RDF::RNode build(const std::string& base_dir, 
                           const VariableRegistry&, 
                           IEventProcessor& processor,
                           const std::string& rel, 
                           const nlohmann::json& all) {
        auto df = this->open_frame(base_dir, rel, processor, origin_);
        df = this->filter_truth(df, truth_);
        df = this->exclude_truth(df, exclude_, all);
        return df;
    }
};

}

#endif
