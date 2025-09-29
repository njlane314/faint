#ifndef ANALYSIS_SAMPLE_SET_H
#define ANALYSIS_SAMPLE_SET_H

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "ROOT/RDataFrame.hxx"

#include "faint/PreSelection.h"
#include "faint/Run.h"
#include "faint/RunCatalog.h"
#include "faint/Sample.h"
#include "faint/SelectionQuery.h"
#include "faint/TruthClassifier.h"
#include "faint/Variables.h"
#include "faint/Weighter.h"

namespace nlohmann {
class json;
}

namespace faint {

class SampleSet {
 public:
  using Map = std::map<SampleKey, Sample>;

  SampleSet(const RunCatalog& runs, VariableRegistry variables,
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
  void snapshot(const SelectionQuery& query, const std::string& out,
                const std::vector<std::string>& cols = {}) const;

  void print_branches();

 private:
  const RunCatalog& runs_;
  VariableRegistry variables_;
  std::string ntuple_dir_;

  std::string beam_;
  std::vector<std::string> periods_;
  bool blind_;

  double total_pot_;
  long total_triggers_;

  Map samples_;
  std::vector<std::unique_ptr<IEventProcessor>> processors_;
  std::unordered_map<SampleKey, const Run*> run_cache_;

  void build();
  void add_run(const Run& rc);
  std::unique_ptr<IEventProcessor> build_pipeline(const nlohmann::json& sample);
  void snapshot_impl(const std::string& filter, const std::string& out,
                     const std::vector<std::string>& cols) const;
};

}  // namespace faint

#endif
