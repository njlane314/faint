#ifndef FAINT_LOG_H
#define FAINT_LOG_H

#include <chrono>
#include <ctime>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <iostream>

namespace faint {
namespace log {

enum class Level { kDebug, kInfo, kWarn, kError };

inline std::mutex& mutex() {
  static std::mutex m;
  return m;
}

inline const char* label(Level level) {
  switch (level) {
    case Level::kDebug:
      return "DEBUG";
    case Level::kInfo:
      return "INFO";
    case Level::kWarn:
      return "WARN";
    case Level::kError:
      return "ERROR";
  }
  return "LOG";
}

template <typename Stream>
void write_header(Stream& os, Level level, const std::string& scope) {
  using clock = std::chrono::system_clock;
  const auto now = clock::now();
  const std::time_t t = clock::to_time_t(now);
  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &t);
#else
  localtime_r(&t, &tm);
#endif
  os << '[' << label(level) << "] " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << " | ";
  if (!scope.empty()) {
    os << scope << ": ";
  }
}

template <typename Stream, typename Arg>
void append(Stream& os, Arg&& arg) {
  os << std::forward<Arg>(arg);
}

template <typename Stream, typename Arg, typename... Args>
void append(Stream& os, Arg&& arg, Args&&... rest) {
  os << std::forward<Arg>(arg);
  if (sizeof...(rest) > 0) os << ' ';
  append(os, std::forward<Args>(rest)...);
}

template <typename... Args>
void log(Level level, const std::string& scope, Args&&... args) {
  std::lock_guard<std::mutex> lock(mutex());
  std::ostream& os = (level == Level::kWarn || level == Level::kError)
                         ? static_cast<std::ostream&>(std::cerr)
                         : static_cast<std::ostream&>(std::clog);
  std::ostringstream buffer;
  append(buffer, std::forward<Args>(args)...);
  write_header(os, level, scope);
  os << buffer.str() << '\n';
}

template <typename... Args>
[[noreturn]] void fatal(const std::string& scope, Args&&... args) {
  std::ostringstream buffer;
  append(buffer, std::forward<Args>(args)...);
  log(Level::kError, scope, buffer.str());
  throw std::runtime_error(buffer.str());
}

template <typename... Args>
void debug(const std::string& scope, Args&&... args) {
  log(Level::kDebug, scope, std::forward<Args>(args)...);
}

template <typename... Args>
void info(const std::string& scope, Args&&... args) {
  log(Level::kInfo, scope, std::forward<Args>(args)...);
}

template <typename... Args>
void warn(const std::string& scope, Args&&... args) {
  log(Level::kWarn, scope, std::forward<Args>(args)...);
}

}  // namespace log
}  // namespace faint

#endif  // FAINT_LOG_H

