#include "faint/Samples.h"

#include <utility>

namespace sample {

SampleKey::SampleKey(std::string value) : value_(std::move(value)) {}

SampleKey::SampleKey(const char* value) : value_(value ? value : "") {}

std::string to_key(SampleVariation var) {
  switch (var) {
    case SampleVariation::kCV:
      return "CV";
    case SampleVariation::kLYAttenuation:
      return "LYAttenuation";
    case SampleVariation::kLYDown:
      return "LYDown";
    case SampleVariation::kLYRayleigh:
      return "LYRayleigh";
    case SampleVariation::kRecomb2:
      return "Recomb2";
    case SampleVariation::kSCE:
      return "SCE";
    case SampleVariation::kWireModX:
      return "WireModX";
    case SampleVariation::kWireModYZ:
      return "WireModYZ";
    case SampleVariation::kWireModAngleXZ:
      return "WireModAngleXZ";
    case SampleVariation::kWireModAngleYZ:
      return "WireModAngleYZ";
    default:
      return "Unknown";
  }
}

}  // namespace sample
