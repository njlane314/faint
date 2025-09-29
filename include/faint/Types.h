#ifndef ANALYSIS_TYPES_H
#define ANALYSIS_TYPES_H

#include <functional>
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

using SampleOrigin = Origin;
using SampleVariation = Variation;

class SampleKey {
  public:
    SampleKey() = default;
    explicit SampleKey(std::string key) : key_(std::move(key)) {}

    const std::string &str() const noexcept { return key_; }

    bool empty() const noexcept { return key_.empty(); }

    friend bool operator==(const SampleKey &lhs, const SampleKey &rhs) noexcept {
        return lhs.key_ == rhs.key_;
    }

    friend bool operator<(const SampleKey &lhs, const SampleKey &rhs) noexcept {
        return lhs.key_ < rhs.key_;
    }

  private:
    std::string key_;
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

namespace std {
template <>
struct hash<faint::SampleKey> {
    std::size_t operator()(const faint::SampleKey &key) const noexcept {
        return std::hash<std::string>{}(key.str());
    }
};
}

#endif
