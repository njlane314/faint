#ifndef ANALYSIS_EVENT_PROCESSOR_H
#define ANALYSIS_EVENT_PROCESSOR_H

#include <memory>

#include "ROOT/RDataFrame.hxx"

namespace faint {
namespace sample {

enum class SampleOrigin : unsigned int;

}  // namespace sample

class EventProcessor {
public:
  virtual ~EventProcessor() = default;

  virtual ROOT::RDF::RNode process(ROOT::RDF::RNode df,
                                   sample::SampleOrigin origin) const = 0;

  void chain_processor(std::unique_ptr<EventProcessor> next);

protected:
  std::unique_ptr<EventProcessor> next_;
};

}

#endif
