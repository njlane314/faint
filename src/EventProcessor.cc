#include "faint/EventProcessor.h"

namespace faint {

void EventProcessor::chain_processor(std::unique_ptr<EventProcessor> next) {
  next_ = std::move(next);
}

}  // namespace faint
