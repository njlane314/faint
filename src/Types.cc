#include "faint/Samples.h"

#include <utility>

namespace faint::sample {

Key::Key(std::string value) : value_(std::move(value)) {}

Key::Key(const char* value) : value_(value ? value : "") {}

std::string to_key(Variation var) {
  switch (var) {
    case Variation::kCV:
      return "CV";
    case Variation::kLYAttenuation:
      return "LYAttenuation";
    case Variation::kLYDown:
      return "LYDown";
    case Variation::kLYRayleigh:
      return "LYRayleigh";
    case Variation::kRecomb2:
      return "Recomb2";
    case Variation::kSCE:
      return "SCE";
    case Variation::kWireModX:
      return "WireModX";
    case Variation::kWireModYZ:
      return "WireModYZ";
    case Variation::kWireModAngleXZ:
      return "WireModAngleXZ";
    case Variation::kWireModAngleYZ:
      return "WireModAngleYZ";
    default:
      return "Unknown";
  }
}

}  // namespace faint::sample
