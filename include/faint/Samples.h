#ifndef FAINT_SAMPLES_H
#define FAINT_SAMPLES_H

#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ROOT/RDataFrame.hxx"

#include "nlohmann/json.hpp"

#include "faint/EventProcessor.h"
#include "faint/PreSelection.h"
#include "faint/Run.h"
#include "faint/Selection.h"
#include "faint/TruthClassifier.h"
#include "faint/Types.h"
#include "faint/Variables.h"
#include "faint/Weighter.h"

namespace faint {

class Sample {
 public:
  Sample(const nlohmann::json& j, const nlohmann::json& all,
         const std::string& base_dir, const VariableRegistry& vars,
         EventProcessor& processor);

  const SampleKey& key() const noexcept { return key_; }
  SampleOrigin origin() const noexcept { return origin_; }
  double pot() const noexcept { return pot_; }
  long triggers() const noexcept { return triggers_; }

  ROOT::RDF::RNode nominal() const { return nominal_node_; }
  const std::map<SampleVariation, ROOT::RDF::RNode>& variations() const {
    return variations_;
  }

 private:
  void validate(const std::string& base_dir) const;

  SampleKey key_;
  SampleOrigin origin_;

  std::string path_;
  std::string truth_;
  std::vector<std::string> exclude_;

  double pot_{0.0};
  long triggers_{0};

  ROOT::RDF::RNode nominal_node_;
  std::map<SampleVariation, ROOT::RDF::RNode> variations_;
  std::map<SampleVariation, std::string> variation_paths_;

  SampleVariation parse_variation(const std::string& s) const;
  ROOT::RDF::RNode build(const std::string& base_dir,
                         [[maybe_unused]] const VariableRegistry& vars,
                         EventProcessor& processor, const std::string& rel,
                         const nlohmann::json& all);
};

class SampleSet {
 public:
  using Map = std::map<SampleKey, Sample>;

  SampleSet(const RunReader& runs, VariableRegistry variables,
            const std::string& beam, std::vector<std::string> periods,
            const std::string& ntuple_dir, bool blind = true);

  Map& frames() noexcept;

  double total_pot() const noexcept;
  long total_triggers() const noexcept;

  const std::string& beam() const noexcept;
  const std::vector<std::string>& periods() const noexcept;

  const Run* run_for(const SampleKey& sk) const;

  void snapshot(const std::string& filter, const std::string& out,
                const std::vector<std::string>& cols = {}) const;
  void snapshot(const Selection& selection, const std::string& out,
                const std::vector<std::string>& cols = {}) const;

  void print_branches();

 private:
  const RunReader& runs_;
  VariableRegistry variables_;
  std::string ntuple_dir_;

  std::string beam_;
  std::vector<std::string> periods_;
  bool blind_;

  double total_pot_;
  long total_triggers_;

  Map samples_;
  std::vector<std::unique_ptr<EventProcessor>> processors_;
  std::unordered_map<SampleKey, const Run*> run_cache_;

  void build();
  void add_run(const Run& rc);
  std::unique_ptr<EventProcessor> build_pipeline(const nlohmann::json& sample);
  void snapshot_impl(const std::string& filter, const std::string& out,
                     const std::vector<std::string>& cols) const;
};

}  // namespace faint

#endif
