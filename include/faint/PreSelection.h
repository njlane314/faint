#ifndef ANALYSIS_NUMUCC_SELECTOR_H
#define ANALYSIS_NUMUCC_SELECTOR_H

#include "ROOT/RVec.hxx"

#include <faint/EventProcessor.h>
#include <faint/Types.h>

namespace faint {

// NuMu charged-current event selection + the useful reconstruction gates
// previously split across Reconstruction and NuMuCCSelectionProcessor.
class PreSelection : public EventProcessor {
public:
  ROOT::RDF::RNode process(ROOT::RDF::RNode df, Origin origin) const override;
};

} // namespace faint

#endif
