#ifndef ANALYSIS_NUMUCC_SELECTOR_H
#define ANALYSIS_NUMUCC_SELECTOR_H

#include "ROOT/RVec.hxx"

#include <rarexsec/EventProcessor.h>
#include <rarexsec/Types.h>

namespace analysis {

// NuMu charged-current event selection + the useful reconstruction gates
// previously split across Reconstruction and NuMuCCSelectionProcessor.
class PreSelection : public EventProcessor {
public:
  ROOT::RDF::RNode process(ROOT::RDF::RNode df, Origin origin) const override;
};

} // namespace analysis

#endif
