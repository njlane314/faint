#ifndef ANALYSIS_TYPES_H
#define ANALYSIS_TYPES_H

#include "faint/Samples.h"

namespace faint {

using SampleKey [[deprecated("Use faint::sample::SampleKey instead")]] =
    sample::SampleKey;
using SampleOrigin [[deprecated("Use faint::sample::SampleOrigin instead")]] =
    sample::SampleOrigin;
using SampleRole [[deprecated("Use faint::sample::SampleRole instead")]] =
    sample::SampleRole;
using SampleVariation [[deprecated("Use faint::sample::SampleVariation instead")]] =
    sample::SampleVariation;

[[deprecated("Use faint::sample::to_key instead")]] inline std::string to_key(
    SampleVariation var) {
  return sample::to_key(var);
}

}  // namespace faint

#endif
