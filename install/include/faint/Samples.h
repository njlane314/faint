#ifndef FAINT_SAMPLES_H
#define FAINT_SAMPLES_H

#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ROOT/RDataFrame.hxx"

#include "nlohmann/json.hpp"

#include <faint/EventProcessor.h>
#include <faint/PreSelection.h>
#include <faint/Run.h>
#include <faint/TruthClassifier.h>
#include <faint/Variables.h>
#include <faint/Weighter.h>

namespace faint::selection {
class Selection;
}

namespace faint::sample {

class Key {
 public:
  Key() = default;
  explicit Key(std::string value);
  Key(const char* value);

  const std::string& str() const noexcept { return value_; }
  const char* c_str() const noexcept { return value_.c_str(); }

  bool empty() const noexcept { return value_.empty(); }

  friend bool operator==(const Key& lhs, const Key& rhs) noexcept {
    return lhs.value_ == rhs.value_;
  }

  friend bool operator!=(const Key& lhs, const Key& rhs) noexcept {
    return !(lhs == rhs);
  }

  friend bool operator<(const Key& lhs, const Key& rhs) noexcept {
    return lhs.value_ < rhs.value_;
  }

 private:
  std::string value_;
};

enum class Origin : unsigned int {
  kUnknown = 0,
  kData,
  kMonteCarlo,
  kExternal,
  kDirt
};

enum class Role { kData, kNominal, kVariation };

enum class Variation : unsigned int {
  kUnknown = 0,
  kCV,
  kLYAttenuation,
  kLYDown,
  kLYRayleigh,
  kRecomb2,
  kSCE,
  kWireModX,
  kWireModYZ,
  kWireModAngleXZ,
  kWireModAngleYZ
};

std::string to_key(Variation var);

}  // namespace faint::sample

namespace faint {

class Sample {
 public:
  Sample(const nlohmann::json& j, const nlohmann::json& all,
         const std::string& base_dir, const VariableRegistry& vars,
         EventProcessor& processor);

  const sample::Key& key() const noexcept { return key_; }
  sample::Origin origin() const noexcept { return origin_; }
  double pot() const noexcept { return pot_; }
  long triggers() const noexcept { return triggers_; }

  ROOT::RDF::RNode nominal() const { return nominal_node_; }
  const std::map<sample::Variation, ROOT::RDF::RNode>& variations() const {
    return variations_;
  }

 private:
  void validate(const std::string& base_dir) const;

  sample::Key key_;
  sample::Origin origin_;

  std::string path_;
  std::string truth_;
  std::vector<std::string> exclude_;

  double pot_{0.0};
  long triggers_{0};

  ROOT::RDF::RNode nominal_node_;
  std::map<sample::Variation, ROOT::RDF::RNode> variations_;
  std::map<sample::Variation, std::string> variation_paths_;

  sample::Variation parse_variation(const std::string& s) const;
  ROOT::RDF::RNode build(const std::string& base_dir,
                         [[maybe_unused]] const VariableRegistry& vars,
                         EventProcessor& processor, const std::string& rel,
                         const nlohmann::json& all);
};

class SampleSet {
 public:
  using Map = std::map<sample::Key, Sample>;

  SampleSet(const RunReader& runs, VariableRegistry variables,
            const std::string& beam, std::vector<std::string> periods,
            const std::string& ntuple_dir);

  Map& frames() noexcept;

  double total_pot() const noexcept;
  long total_triggers() const noexcept;

  const std::string& beam() const noexcept;
  const std::vector<std::string>& periods() const noexcept;

  const Run* run_for(const sample::Key& sk) const;

  void snapshot(const std::string& filter, const std::string& out,
                const std::vector<std::string>& cols = {}) const;
  void snapshot(const selection::Selection& selection, const std::string& out,
                const std::vector<std::string>& cols = {}) const;

  void print_branches();

 private:
  const RunReader& runs_;
  VariableRegistry variables_;
  std::string ntuple_dir_;

  std::string beam_;
  std::vector<std::string> periods_;
  double total_pot_;
  long total_triggers_;

  Map samples_;
  std::vector<std::unique_ptr<EventProcessor>> processors_;
  std::unordered_map<sample::Key, const Run*> run_cache_;

  void build();
  void add_run(const Run& rc);
  std::unique_ptr<EventProcessor> build_pipeline(const nlohmann::json& sample);
  void snapshot_impl(const std::string& filter, const std::string& out,
                     const std::vector<std::string>& cols) const;
};

}  // namespace faint

namespace std {

template <>
struct hash<faint::sample::Key> {
  std::size_t operator()(const faint::sample::Key& key) const noexcept {
    return std::hash<std::string>{}(key.str());
  }
};

}  // namespace std

#endif
