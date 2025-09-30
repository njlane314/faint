#ifndef ANALYSIS_TYPES_H
#define ANALYSIS_TYPES_H

#include <functional>
#include <string>
#include <utility>

namespace faint {

class SampleKey {
 public:
  SampleKey() = default;
  explicit SampleKey(std::string value);
  SampleKey(const char* value);

  const std::string& str() const noexcept { return value_; }
  const char* c_str() const noexcept { return value_.c_str(); }

  bool empty() const noexcept { return value_.empty(); }

  friend bool operator==(const SampleKey& lhs, const SampleKey& rhs) noexcept {
    return lhs.value_ == rhs.value_;
  }

  friend bool operator!=(const SampleKey& lhs, const SampleKey& rhs) noexcept {
    return !(lhs == rhs);
  }

  friend bool operator<(const SampleKey& lhs, const SampleKey& rhs) noexcept {
    return lhs.value_ < rhs.value_;
  }

 private:
  std::string value_;
};

enum class SampleOrigin : unsigned int {
  kUnknown = 0,
  kData,
  kMonteCarlo,
  kExternal,
  kDirt
};

enum class SampleRole { kData, kNominal, kVariation };

enum class SampleVariation : unsigned int {
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

std::string to_key(SampleVariation var);

}  // namespace faint

namespace std {

template <>
struct hash<faint::SampleKey> {
  std::size_t operator()(const faint::SampleKey& key) const noexcept {
    return std::hash<std::string>{}(key.str());
  }
};

}  // namespace std

#endif
