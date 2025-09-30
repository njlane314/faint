#ifndef ANALYSIS_EVENT_PROCESSOR_H
#define ANALYSIS_EVENT_PROCESSOR_H

#include <memory>

#include "ROOT/RDataFrame.hxx"

namespace sample {
enum class SampleOrigin : unsigned int;
}

namespace faint {

using sample::SampleOrigin;

class EventProcessor {
public:
  virtual ~EventProcessor() = default;

  virtual ROOT::RDF::RNode process(ROOT::RDF::RNode df,
                                   SampleOrigin origin) const = 0;

  void chain_processor(std::unique_ptr<EventProcessor> next);

protected:
  std::unique_ptr<EventProcessor> next_;
};

} 

#endif
