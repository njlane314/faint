#ifndef FAINT_LOGGER_H
#define FAINT_LOGGER_H

#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace faint {
namespace log {

enum class Level {
    kDebug = 0,
    kInfo,
    kWarn,
    kError,
    kFatal
};

namespace detail {

inline std::mutex &logger_mutex() {
    static std::mutex mutex;
    return mutex;
}

inline Level &current_level() {
    static Level level = Level::kInfo;
    return level;
}

inline std::string level_to_string(Level level) {
    switch (level) {
        case Level::kDebug: return "DEBUG";
        case Level::kInfo: return "INFO";
        case Level::kWarn: return "WARN";
        case Level::kError: return "ERROR";
        case Level::kFatal: return "FATAL";
    }
    return "UNKNOWN";
}

inline std::string timestamp() {
    using clock = std::chrono::system_clock;
    const auto now = clock::now();
    const auto time = clock::to_time_t(now);

    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &time);
#else
    localtime_r(&time, &tm);
#endif

    std::ostringstream ss;
    ss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

template <typename... Args>
inline std::string build_message(Args &&...args) {
    std::ostringstream ss;
    (ss << ... << std::forward<Args>(args));
    return ss.str();
}

inline std::ostream &select_stream(Level level) {
    if (level == Level::kWarn || level == Level::kError || level == Level::kFatal) {
        return std::cerr;
    }
    return std::cout;
}

template <typename... Args>
inline void write(Level level, Args &&...args) {
    if (level < current_level()) return;

    const std::string message = build_message(std::forward<Args>(args)...);

    std::lock_guard<std::mutex> lock(logger_mutex());
    auto &stream = select_stream(level);
    stream << "[" << timestamp() << "] [" << level_to_string(level) << "] " << message << std::endl;

    if (level == Level::kFatal) {
        throw std::runtime_error(message);
    }
}

}  // namespace detail

inline void set_level(Level level) { detail::current_level() = level; }

inline Level level() { return detail::current_level(); }

template <typename... Args>
inline void debug(Args &&...args) {
    detail::write(Level::kDebug, std::forward<Args>(args)...);
}

template <typename... Args>
inline void info(Args &&...args) {
    detail::write(Level::kInfo, std::forward<Args>(args)...);
}

template <typename... Args>
inline void warn(Args &&...args) {
    detail::write(Level::kWarn, std::forward<Args>(args)...);
}

template <typename... Args>
inline void error(Args &&...args) {
    detail::write(Level::kError, std::forward<Args>(args)...);
}

template <typename... Args>
inline void fatal(Args &&...args) {
    detail::write(Level::kFatal, std::forward<Args>(args)...);
}

}  // namespace log
}  // namespace faint

#endif  // FAINT_LOGGER_H

