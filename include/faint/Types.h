#ifndef ANALYSIS_TYPES_H
#define ANALYSIS_TYPES_H

#include <string>

namespace analysis {

enum class Origin : unsigned int { kUnknown = 0, kData, kMonteCarlo, kExternal, kDirt };
enum class Role { kData, kNominal, kVariation };

enum class Variation : unsigned int {
    kUnknown = 0,
    kCV,
    kLYAttenuation,
    kLYDown,
    kLYRayleigh,
    kRecomb2,
    kSCE,
    kWireModX,
    kWireModYZ,
    kWireModAngleXZ,
    kWireModAngleYZ
};

inline std::string to_key(Variation var) {
    switch (var) {
        case Variation::kCV:              return "CV";
        case Variation::kLYAttenuation:   return "LYAttenuation";
        case Variation::kLYDown:          return "LYDown";
        case Variation::kLYRayleigh:      return "LYRayleigh";
        case Variation::kRecomb2:         return "Recomb2";
        case Variation::kSCE:             return "SCE";
        case Variation::kWireModX:        return "WireModX";
        case Variation::kWireModYZ:       return "WireModYZ";
        case Variation::kWireModAngleXZ:  return "WireModAngleXZ";
        case Variation::kWireModAngleYZ:  return "WireModAngleYZ";
        default:                          return "Unknown";
    }
}

} 

#endif
