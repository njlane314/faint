#ifndef ANALYSIS_MUON_SELECTOR_H
#define ANALYSIS_MUON_SELECTOR_H

#include <cmath>

#include "ROOT/RVec.hxx"

#include <rarexsec/EventProcessor.h>

namespace analysis {

class MuonSelector : public EventProcessor {
public:
  ROOT::RDF::RNode process(ROOT::RDF::RNode df, Origin origin) const override;

private:
  ROOT::RDF::RNode build_mask(ROOT::RDF::RNode df) const;
  ROOT::RDF::RNode extract_features(ROOT::RDF::RNode df) const;
};

} // namespace analysis

#endif
