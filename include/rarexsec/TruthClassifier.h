#ifndef ANALYSIS_TRUTH_CLASSIFIER_H
#define ANALYSIS_TRUTH_CLASSIFIER_H

#include <cmath>
#include <iostream>
#include <map>
#include <mutex>

#include <rarexsec/EventProcessor.h>

namespace analysis {

class TruthClassifier : public EventProcessor {
public:
  ROOT::RDF::RNode process(ROOT::RDF::RNode df, Origin origin) const override;

private:
  ROOT::RDF::RNode processNonMc(ROOT::RDF::RNode df, Origin origin) const;
  ROOT::RDF::RNode defineCounts(ROOT::RDF::RNode df) const;
  ROOT::RDF::RNode assignInclusiveChannels(ROOT::RDF::RNode df) const;
  ROOT::RDF::RNode assignExclusiveChannels(ROOT::RDF::RNode df) const;
  ROOT::RDF::RNode assignChannelDefinitions(ROOT::RDF::RNode df) const;
};

// ---------- Backward-compat (remove when migrated) ----------
using TruthChannelProcessor = TruthClassifier;

} // namespace analysis

#endif
