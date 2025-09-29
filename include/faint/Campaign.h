#ifndef FAINT_CAMPAIGN_H
#define FAINT_CAMPAIGN_H

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <utility>

#include "ROOT/RDataFrame.hxx"
#include "TSystem.h"

#include <faint/Types.h>
#include <faint/Variables.h>
#include <faint/RunCatalog.h>
#include <faint/Sample.h>
#include <faint/SampleSet.h>

namespace faint {
namespace campaign {

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

std::string run_config_path();

std::string ntuple_directory();

struct Options {
    std::string beam;
    std::vector<std::string> periods;
    std::string ntuple_dir;
    bool blind{true};
};

class Campaign {
public:
    static Campaign open(const std::string& run_config_json, Options opt, Variables vars = Variables{});

    std::vector<std::string> sample_keys(
        SampleOrigin origin_filter = SampleOrigin::kUnknown) const;

    ROOT::RDF::RNode df(std::string_view sample_key,
                        SampleVariation v = SampleVariation::kCV) const;

    ROOT::RDF::RNode final(std::string_view key,
                           SampleVariation v = SampleVariation::kCV) const;

    ROOT::RDF::RNode quality(std::string_view key,
                             SampleVariation v = SampleVariation::kCV) const;

    void snapshot_where(const std::string& filter, const std::string& out_file, const std::vector<std::string>& columns = {}) const;

    void snapshot_final(const std::string& out_file, const std::vector<std::string>& columns = {}) const;

    double pot() const noexcept;

    long triggers() const noexcept;

    const std::string& beam() const noexcept;

    const std::vector<std::string>& periods() const noexcept;

    const SampleSet& samples() const;

    const RunCatalog& runs() const;

private:
    RunCatalog runs_;
    Variables vars_;
    Options opt_;
    std::unique_ptr<SampleSet> set_;

    const SampleSet& set() const;

    const Sample* find_sample(std::string_view key) const;
};

} // namespace campaign
} // namespace faint

#endif
