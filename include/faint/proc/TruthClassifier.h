#ifndef ANALYSIS_TRUTH_CLASSIFIER_H
#define ANALYSIS_TRUTH_CLASSIFIER_H

#include <cmath>
#include <iostream>
#include <map>
#include <mutex>

#include <faint/proc/EventProcessor.h>

namespace faint {

class TruthClassifier : public EventProcessor {
public:
  ROOT::RDF::RNode process(ROOT::RDF::RNode df,
                           sample::Origin origin) const override;

private:
  ROOT::RDF::RNode process_non_mc(ROOT::RDF::RNode df,
                                sample::Origin origin) const;
  ROOT::RDF::RNode define_counts(ROOT::RDF::RNode df) const;
  ROOT::RDF::RNode assign_inclusive_channels(ROOT::RDF::RNode df) const;
  ROOT::RDF::RNode assign_exclusive_channels(ROOT::RDF::RNode df) const;
  ROOT::RDF::RNode assign_channel_definitions(ROOT::RDF::RNode df) const;
};

} // namespace faint

#endif
