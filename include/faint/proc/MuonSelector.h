#ifndef ANALYSIS_MUON_SELECTOR_H
#define ANALYSIS_MUON_SELECTOR_H

#include <cmath>

#include "ROOT/RVec.hxx"

#include <faint/proc/EventProcessor.h>

namespace faint {

class MuonSelector : public EventProcessor {
public:
  ROOT::RDF::RNode process(ROOT::RDF::RNode df,
                           SampleOrigin origin) const override;

private:
  ROOT::RDF::RNode build_mask(ROOT::RDF::RNode df) const;
  ROOT::RDF::RNode extract_features(ROOT::RDF::RNode df) const;
};

} // namespace faint

#endif
