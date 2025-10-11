#pragma once

namespace rarexsec {
namespace fiducial {

inline constexpr float min_x = 5.f;
inline constexpr float max_x = 251.f;
inline constexpr float min_y = -110.f;
inline constexpr float max_y = 110.f;
inline constexpr float min_z = 20.f;
inline constexpr float max_z = 986.f;

inline constexpr float reco_gap_min_z = 675.f;
inline constexpr float reco_gap_max_z = 775.f;

namespace detail {

template <typename T>
inline bool is_within(const T& value, float low, float high) {
    return value > low && value < high;
}

template <typename X, typename Y, typename Z>
inline bool is_in_active_volume(const X& x, const Y& y, const Z& z) {
    return is_within(x, min_x, max_x) &&
           is_within(y, min_y, max_y) &&
           is_within(z, min_z, max_z);
}

}

template <typename X, typename Y, typename Z>
inline bool is_in_truth_volume(const X& x, const Y& y, const Z& z) {
    return detail::is_in_active_volume(x, y, z);
}

template <typename X, typename Y, typename Z>
inline bool is_in_reco_volume(const X& x, const Y& y, const Z& z) {
    return detail::is_in_active_volume(x, y, z) && (z < reco_gap_min_z || z > reco_gap_max_z);
}

}
}