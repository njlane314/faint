#pragma once

#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "rarexsec/Hub.h"

namespace rarexsec {
struct Env {
  std::string cfg, beamline;
  std::vector<std::string> periods;
  static Env from_env() {
    auto get_env = [](const char* key) {
      const char* value = std::getenv(key);
      return value ? std::string(value) : std::string();
    };

    Env env;
    env.cfg = get_env("RAREXSEC_CFG");
    if (env.cfg.empty()) {
      throw std::runtime_error("RAREXSEC_CFG missing");
    }
    env.beamline = get_env("RAREXSEC_BEAMLINE");
    if (env.beamline.empty()) {
      throw std::runtime_error("RAREXSEC_BEAMLINE missing");
    }
    auto periods = get_env("RAREXSEC_PERIODS");
    if (periods.empty()) {
      throw std::runtime_error("RAREXSEC_PERIODS missing");
    }
    std::replace(periods.begin(), periods.end(), ',', ' ');
    std::stringstream ss(periods);
    std::string token;
    while (ss >> token) {
      env.periods.push_back(token);
    }
    return env;
  }
  Hub make_hub() const { return Hub(cfg); }
};
} // namespace rarexsec
