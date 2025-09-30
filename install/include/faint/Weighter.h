#ifndef ANALYSIS_WEIGHTER_H
#define ANALYSIS_WEIGHTER_H

#include <cmath>
#include <nlohmann/json.hpp>

#include <faint/EventProcessor.h>

namespace faint {

class Weighter : public EventProcessor {
public:
  Weighter(const nlohmann::json& cfg, double total_run_pot,
           long total_run_triggers);

  ROOT::RDF::RNode process(ROOT::RDF::RNode df,
                           SampleOrigin origin) const override;

private:
  double sample_pot_;
  long   sample_triggers_;
  double total_run_pot_;
  long   total_run_triggers_;
};
  
} // namespace faint

#endif
