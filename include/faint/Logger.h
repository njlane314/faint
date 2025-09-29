#ifndef FAINT_LOGGER_H
#define FAINT_LOGGER_H

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace faint {
namespace log {
namespace detail {
inline void append_to_stream(std::ostream&) {}

template <typename T, typename... Rest>
void append_to_stream(std::ostream& os, T&& value, Rest&&... rest) {
  os << std::forward<T>(value);
  if constexpr (sizeof...(rest) > 0) {
    os << ' ';
    append_to_stream(os, std::forward<Rest>(rest)...);
  }
}

template <typename... Args>
std::string build_message(Args&&... args) {
  std::ostringstream ss;
  append_to_stream(ss, std::forward<Args>(args)...);
  return ss.str();
}
}  // namespace detail

template <typename... Args>
void debug(Args&&... args) {
  std::clog << "[DEBUG] " << detail::build_message(std::forward<Args>(args)...) << std::endl;
}

template <typename... Args>
void info(Args&&... args) {
  std::clog << "[INFO] " << detail::build_message(std::forward<Args>(args)...) << std::endl;
}

template <typename... Args>
void warn(Args&&... args) {
  std::clog << "[WARN] " << detail::build_message(std::forward<Args>(args)...) << std::endl;
}

template <typename... Args>
[[noreturn]] void fatal(Args&&... args) {
  throw std::runtime_error(detail::build_message(std::forward<Args>(args)...));
}

}  // namespace log
}  // namespace faint

#endif  // FAINT_LOGGER_H
