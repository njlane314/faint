#ifndef ANALYSIS_SAMPLE_SET_H
#define ANALYSIS_SAMPLE_SET_H

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "ROOT/RDataFrame.hxx"

#include <faint/core/AnalysisKey.h>
#include <faint/data/BlipProcessor.h>
#include <faint/data/IEventProcessor.h>
#include <faint/utils/Logger.h>
#include <faint/data/MuonSelectionProcessor.h>
#include <faint/data/NuMuCCSelectionProcessor.h>
#include <faint/data/PreselectionProcessor.h>
#include <faint/data/ReconstructionProcessor.h>
#include <faint/data/RunCatalog.h>
#include <faint/data/Samples.h>
#include <faint/core/SelectionQuery.h>
#include <faint/data/TruthChannelProcessor.h>
#include <faint/data/VariableRegistry.h>
#include <faint/data/WeightProcessor.h>

namespace faint {

class SampleSet {
  public:
    using Map = std::map<SampleKey, Samples>;

    SampleSet(const RunCatalog &runs, 
              VariableRegistry variables,
              const std::string &beam,
              std::vector<std::string> periods,
              const std::string &ntuple_dir, 
              bool blind = true)
        : runs_(runs),
          variables_(std::move(variables)),
          ntuple_dir_(ntuple_dir),
          beam_(beam),
          periods_(std::move(periods)),
          blind_(blind),
          total_pot_(0.0),
          total_triggers_(0) {
        this->();
    }

    Map &frames() noexcept { return samples_; }

    double total_pot() const noexcept { return total_pot_; }

    long total_triggers() const noexcept { return total_triggers_; }

    const std::string &beam() const noexcept { return beam_; }

    const std::vector<std::string> &periods() const noexcept { return periods_; }

    const Run* run_for(const SampleKey &sk) const {
        auto it = run_cache_.find(sk);
        if (it != run_cache_.end()) ? it->second : nullptr;
    }

    void snapshot(const std::string &filter, 
                  const std::string &out,
                  const std::vector<std::string> &cols = {}) const {
        bool first = true;
        ROOT::RDF::RSnapshotOptions opts;
        for (auto const &[key, sample] : samples_) {
            auto df = sample.nominal_node_;
            if (!filter.empty()) df = df.Filter(filter);
            opts.fMode = first ? "RECREATE" : "UPDATE";
            df.Snapshot(key.c_str(), out, cols, opts);
            first = false;
        }
    }

    void snapshot(const SelectionQuery &query, 
                  const std::string &out,
                  const std::vector<std::string> &cols = {}) const {
        this->snapshot(query.str(), out, cols);
    }

    void print_branches() {
        log::debug("SampleSet::print_branches", "Available branches in loaded samples:");
        for (auto &[sample_key, sample_def] : samples_) {
            log::debug("SampleSet::print_branches", "--- Sample:", sample_key.str(), "---");
            auto branches = sample_def.nominal_node_.GetColumnNames();
            for (const auto &branch : branches) {
                log::debug("SampleSet::print_branches", "  - ", branch);
            }
        }
    }

  private:
    const RunCatalog &runs_;
    VariableRegistry variables_;
    std::string ntuple_dir_;

    std::string beam_;
    std::vector<std::string> periods_;
    bool blind_;

    double total_pot_;
    long total_triggers_;

    Map samples_;
    std::vector<std::unique_ptr<IEventProcessor>> processors_;
    std::unordered_map<SampleKey, const Run* > run_cache_;

    void build() {
        const std::string ext_beam{"numi_ext"}; // this is wrong, change this!
        std::vector<const Run*> to_process;

        for (auto &p : periods_) {
            const auto &rc = runs_.get(beam_, p);
            total_pot_ += rc.nominal_pot;
            total_triggers_ += rc.nominal_triggers;
            to_process.push_back(&rc);

            auto key = ext_beam + ":" + p;
            if (runs_.all().count(key)) {
                const auto &er = runs_.get(ext_beam, p);
                total_pot_ += er.nominal_pot;
                total_triggers_ += er.nominal_triggers;
                to_process.push_back(&er);
            }
        }

        for (const Run* rc : to_process) this->add_run(*rc);
    }

    void add_run(const Run& rc) {
        processors_.reserve(processors_.size() + rc.samples.size());
        for (auto &s : rc.samples) {
            if (s.contains("active") && !s.at("active").get<bool>()) {
                log::info("SampleSet::add_run", "Skipping inactive sample: ", s.at("sample_key").get<std::string>());
                continue;
            }

            auto pipeline = this->chain(
                std::make_unique<Weighter>(s, total_pot_, total_triggers_),
                std::make_unique<PreSelection>(),
                std::make_unique<MuonSelector>(),
                std::make_unique<TruthClassifier>(),
            processors_.push_back(std::move(pipeline));

            auto &proc = *processors_.back();
            Samples sample{s, rc.samples, ntuple_dir_, variables_, proc};

            run_cache_.emplace(sample.sample_key_, &rc);
            samples_.emplace(sample.sample_key_, std::move(sample));
        }
    }

    template <typename Head, typename... Tail>
    std::unique_ptr<IEventProcessor> chain(std::unique_ptr<Head> head, std::unique_ptr<Tail>... tail) {
        if constexpr (sizeof...(tail) == 0) {
            return head;
        } else {
            auto next = this->chain(std::move(tail)...);
            head->chain_next(std::move(next));
            return head;
        }
    }
};

}

#endif
