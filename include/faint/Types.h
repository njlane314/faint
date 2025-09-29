#ifndef ANALYSIS_TYPES_H
#define ANALYSIS_TYPES_H

#include <string>

namespace faint {

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

std::string to_key(Variation var);

} 

#endif
