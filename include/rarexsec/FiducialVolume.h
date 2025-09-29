#ifndef ANALYSIS_FIDUCIALVOLUME_H
#define ANALYSIS_FIDUCIALVOLUME_H

namespace analysis {

namespace fiducial {

inline constexpr float kMinX = 5.f;
inline constexpr float kMaxX = 251.f;
inline constexpr float kMinY = -110.f;
inline constexpr float kMaxY = 110.f;
inline constexpr float kMinZ = 20.f;
inline constexpr float kMaxZ = 986.f;
inline constexpr float kRecoGapMinZ = 675.f;
inline constexpr float kRecoGapMaxZ = 775.f;

namespace detail {

template <typename T>
inline bool is_within(const T &value, float low, float high) {
  return value > low && value < high;
}

template <typename X, typename Y, typename Z>
inline bool is_in_active_volume(const X &x, const Y &y, const Z &z) {
  return is_within(x, kMinX, kMaxX) && is_within(y, kMinY, kMaxY) &&
         is_within(z, kMinZ, kMaxZ);
}

} // namespace detail

template <typename X, typename Y, typename Z>
inline bool is_in_truth_volume(const X &x, const Y &y, const Z &z) {
  return detail::is_in_active_volume(x, y, z);
}

template <typename X, typename Y, typename Z>
inline bool is_in_reco_volume(const X &x, const Y &y, const Z &z) {
  return detail::is_in_active_volume(x, y, z) &&
         (z < kRecoGapMinZ || z > kRecoGapMaxZ);
}

} // namespace fiducial

} // namespace analysis

#endif
