#ifndef FAINT_DATASET_H
#define FAINT_DATASET_H

#include <algorithm>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "TSystem.h"

#include <faint/Types.h>
#include <faint/Variables.h>
#include <faint/Run.h>
#include <faint/Samples.h>
#include <faint/Selections.h>

namespace faint {
namespace dataset {

namespace sel {
inline constexpr const char* Pre = selection::column::kPassPre;
inline constexpr const char* Flash = selection::column::kPassFlash;
inline constexpr const char* FV = selection::column::kPassFiducial;
inline constexpr const char* Muon = selection::column::kPassMuon;
inline constexpr const char* Topo = selection::column::kPassTopology;
inline constexpr const char* Final = selection::column::kPassFinal;
inline constexpr const char* Quality = selection::column::kQualityEvent;
}

namespace col {
inline constexpr const char* Weight = "nominal_event_weight";
}

std::string run_config_path();

std::string ntuple_directory();

std::string ntuple_directory(const std::string& run_config_json);

struct Options {
    std::string beam;
    std::vector<std::string> periods;
    std::string ntuple_dir;
    bool blind{true};
};

class Dataset {
public:
    struct Entry {
        SampleOrigin origin{SampleOrigin::kUnknown};
        SampleRole role{SampleRole::kNominal};

        Entry() = default;

        Entry(SampleOrigin origin, SampleRole role, ROOT::RDF::RNode node)
            : origin(origin), role(role), dataframe_(std::move(node)) {}

        ROOT::RDF::RNode dataframe() const {
            if (!dataframe_) {
                throw std::runtime_error("Dataset entry does not hold a dataframe");
            }
            return *dataframe_;
        }

    private:
        mutable std::optional<ROOT::RDF::RNode> dataframe_;
    };

    struct Variations {
        Entry nominal;
        std::unordered_map<SampleVariation, Entry> variations;
    };

    using Map = std::unordered_map<SampleKey, Variations>;

    static Dataset open(const std::string& run_config_json, Options opt,
                        Variables vars = Variables{});

    std::vector<std::string> sample_keys(
        SampleOrigin origin_filter = SampleOrigin::kUnknown) const;

    ROOT::RDF::RNode df(std::string_view sample_key,
                        SampleVariation v = SampleVariation::kCV) const;

    ROOT::RDF::RNode final(std::string_view key,
                           SampleVariation v = SampleVariation::kCV) const;

    ROOT::RDF::RNode quality(std::string_view key,
                             SampleVariation v = SampleVariation::kCV) const;

    void snapshot_where(const std::string& filter, const std::string& out_file,
                        const std::vector<std::string>& columns = {}) const;

    void snapshot_final(const std::string& out_file,
                        const std::vector<std::string>& columns = {}) const;

    double pot() const noexcept;

    long triggers() const noexcept;

    const std::string& beam() const noexcept;

    const std::vector<std::string>& periods() const noexcept;

    const SampleSet& samples() const;

    const RunReader& runs() const;

    const Map& datasets() const noexcept { return datasets_; }

private:
    RunReader runs_;
    Variables vars_;
    Options opt_;
    std::unique_ptr<SampleSet> set_;
    Map datasets_;

    const SampleSet& set() const;

    void build_dataset_cache();

    const Variations* find_dataset(std::string_view key) const;

    static Entry make_entry(const Sample& sample, SampleVariation variation,
                            ROOT::RDF::RNode node);
};

} // namespace dataset
} // namespace faint

#endif
