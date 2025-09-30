#ifndef ANALYSIS_NUMUCC_SELECTOR_H
#define ANALYSIS_NUMUCC_SELECTOR_H

#include "ROOT/RVec.hxx"

#include <faint/EventProcessor.h>
#include <faint/Types.h>

namespace faint {

class PreSelection : public EventProcessor {
public:
  ROOT::RDF::RNode process(ROOT::RDF::RNode df,
                           SampleOrigin origin) const override;
};

} // namespace faint

#endif
