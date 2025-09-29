#ifndef SAMPLE_DEFINITION_H
#define SAMPLE_DEFINITION_H

#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "nlohmann/json.hpp"

#include "faint/EventProcessor.h"
#include "faint/Types.h"
#include "faint/Variables.h"

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

}  // namespace faint

#endif
