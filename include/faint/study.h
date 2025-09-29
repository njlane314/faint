#ifndef FAINT_STUDY_H
#define FAINT_STUDY_H

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "TSystem.h"

#include <faint/Types.h>
#include <faint/Variables.h>
#include <faint/data/RunCatalog.h>
#include <faint/data/Sample.h>
#include <faint/data/SampleSet.h>

namespace analysis {
namespace study {

namespace sel {
inline constexpr const char* Pre = "pass_pre";
inline constexpr const char* Flash = "pass_flash";
inline constexpr const char* FV = "pass_fv";
inline constexpr const char* Muon = "pass_mu";
inline constexpr const char* Topo = "pass_topo";
inline constexpr const char* Final = "pass_final";
inline constexpr const char* Quality = "quality_event";
}

namespace col {
inline constexpr const char* Weight = "nominal_event_weight";
}

inline std::string run_config_path() {
    const char* env = gSystem->Getenv("FAINT_RUN_CONFIG");
    if (env) return env;
    std::string cwd = gSystem->pwd();
    return cwd + "/data/samples.json";
}

inline std::string ntuple_directory() {
    const char* env = gSystem->Getenv("FAINT_NTUPLES");
    if (!env) throw std::runtime_error("Set FAINT_NTUPLES to the directory containing faint ntuples");
    return env;
}

struct Options {
    std::string beam;
    std::vector<std::string> periods;
    std::string ntuple_dir;
    bool blind{true};
};

class Study {
public:
    static Study open(const std::string& run_config_json, Options opt, Variables vars = Variables{}) {
        Study s;
        s.runs_ = RunCatalog::fromFile(run_config_json);
        s.vars_ = std::move(vars);
        s.opt_ = std::move(opt);
        s.set_ = std::make_unique<SampleSet>(s.runs_, s.vars_, s.opt_.beam, s.opt_.periods, s.opt_.ntuple_dir, s.opt_.blind);
        return s;
    }

    std::vector<std::string> sample_keys(Origin origin_filter = Origin::kUnknown) const {
        std::vector<std::string> out;
        for (auto const& kv : set().frames()) {
            const auto& key = kv.first;
            const auto& s = kv.second;
            if (origin_filter == Origin::kUnknown || s.origin_ == origin_filter) out.push_back(key.str());
        }
        std::sort(out.begin(), out.end());
        return out;
    }

    ROOT::RDF::RNode df(std::string_view sample_key, Variation v = Variation::kCV) const {
        const Sample* s = find_sample(sample_key);
        if (!s) throw std::runtime_error(std::string("Sample not found: ") + std::string(sample_key));
        if (v == Variation::kCV) return s->node_;
        auto it = s->variations_.find(v);
        return (it != s->variations_.end()) ? it->second : s->node_;
    }

    ROOT::RDF::RNode final(std::string_view key, Variation v = Variation::kCV) const {
        return df(key, v).Filter(sel::Final);
    }

    ROOT::RDF::RNode quality(std::string_view key, Variation v = Variation::kCV) const {
        return df(key, v).Filter(sel::Quality);
    }

    void snapshot_where(const std::string& filter, const std::string& out_file, const std::vector<std::string>& columns = {}) const {
        set().snapshot(filter, out_file, columns);
    }

    void snapshot_final(const std::string& out_file, const std::vector<std::string>& columns = {}) const {
        snapshot_where(sel::Final, out_file, columns);
    }

    double pot() const noexcept { return set().total_pot(); }

    long triggers() const noexcept { return set().total_triggers(); }

    const std::string& beam() const noexcept { return opt_.beam; }

    const std::vector<std::string>& periods() const noexcept { return opt_.periods; }

    const SampleSet& samples() const { return set(); }

    const RunCatalog& runs() const { return runs_; }

private:
    RunCatalog runs_;
    Variables vars_;
    Options opt_;
    std::unique_ptr<SampleSet> set_;

    const SampleSet& set() const { return *set_; }

    const Sample* find_sample(std::string_view key) const {
        for (auto const& kv : set().frames()) {
            if (kv.first.str() == key) return &kv.second;
        }
        return nullptr;
    }
};

}
}

#endif
