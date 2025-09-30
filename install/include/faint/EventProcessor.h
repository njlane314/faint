#ifndef ANALYSIS_EVENT_PROCESSOR_H
#define ANALYSIS_EVENT_PROCESSOR_H

#include <memory>

#include "ROOT/RDataFrame.hxx"

namespace faint::sample {
enum class Origin : unsigned int;
}

namespace faint {

class EventProcessor {
public:
  virtual ~EventProcessor() = default;

  virtual ROOT::RDF::RNode process(ROOT::RDF::RNode df,
                                   sample::Origin origin) const = 0;

  void chain_processor(std::unique_ptr<EventProcessor> next);

protected:
  std::unique_ptr<EventProcessor> next_;
};

} 

#endif
