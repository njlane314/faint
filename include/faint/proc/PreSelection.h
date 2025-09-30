#ifndef ANALYSIS_NUMUCC_SELECTOR_H
#define ANALYSIS_NUMUCC_SELECTOR_H

#include "ROOT/RVec.hxx"

#include <faint/proc/EventProcessor.h>

namespace faint {

class PreSelection : public EventProcessor {
public:
  ROOT::RDF::RNode process(ROOT::RDF::RNode df,
                           sample::SampleOrigin origin) const override;
};

} // namespace faint

#endif
