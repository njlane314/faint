#ifndef ANALYSIS_EVENT_PROCESSOR_H
#define ANALYSIS_EVENT_PROCESSOR_H

#include <memory>

#include "ROOT/RDataFrame.hxx"

#include <faint/Types.h>

namespace faint {

class EventProcessor {
public:
  virtual ~EventProcessor() = default;

  virtual ROOT::RDF::RNode process(ROOT::RDF::RNode df, Origin origin) const = 0;

  void chain_processor(std::unique_ptr<EventProcessor> next) { next_ = std::move(next); }

protected:
  std::unique_ptr<EventProcessor> next_;
};

} 

#endif
