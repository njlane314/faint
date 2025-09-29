#ifndef SAMPLE_DEFINITION_H
#define SAMPLE_DEFINITION_H

#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "nlohmann/json.hpp"

#include <faint/data/IEventProcessor.h>
#include <faint/data/SampleTypes.h>
#include <faint/data/VariableRegistry.h>

namespace faint {

class Sample {
 public:
  Sample(const nlohmann::json& j, const nlohmann::json& all,
         const std::string& base_dir, const VariableRegistry& vars,
         IEventProcessor& processor);

  void validate(const std::string& base_dir) const;

  SampleKey key_;
  SampleOrigin origin_;

  std::string path_;
  std::string truth_;
  std::vector<std::string> exclude_;

  double pot_{0.0};
  long triggers_{0};

  ROOT::RDF::RNode node_;
  std::map<SampleVariation, ROOT::RDF::RNode> variations_;

 private:
  std::map<SampleVariation, std::string> variation_paths_;

  SampleVariation parse_variation(const std::string& s) const;
  ROOT::RDF::RNode build(const std::string& base_dir, const VariableRegistry& vars,
                         IEventProcessor& processor, const std::string& rel,
                         const nlohmann::json& all);
};

}  // namespace faint

#endif
