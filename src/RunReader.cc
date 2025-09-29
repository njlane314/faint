#include "faint/RunReader.h"

#include <fstream>
#include <stdexcept>
#include <utility>

#include "faint/utils/Logger.h"

namespace faint {

void RunReader::add(Run rc) {
  auto key = rc.label();
  if (configs_.count(key))
    throw std::runtime_error("Duplicate Run label: " + key);
  configs_.emplace(std::move(key), std::move(rc));
}

const Run& RunReader::get(const std::string& beam,
                          const std::string& period) const {
  auto key = beam + ":" + period;
  auto it = configs_.find(key);
  if (it == configs_.end())
    throw std::out_of_range("Run not found: " + key);
  return it->second;
}

const std::map<std::string, Run>& RunReader::all() const noexcept {
  return configs_;
}

RunReader RunReader::read(const nlohmann::json& data) {
  RunReader out;
  const std::string top =
      data.contains("run_configurations") ? "run_configurations" : "beamlines";
  for (auto const& [beam, runs] : data.at(top).items()) {
    for (auto const& [period, details] : runs.items()) {
      Run rc(details, beam, period);
      rc.validate();
      out.add(std::move(rc));
    }
  }
  return out;
}

RunReader RunReader::read(const std::string& path) {
  std::ifstream f(path);
  if (!f.is_open())
    log::fatal("RunReader::read", "Could not open config file:", path);
  try {
    nlohmann::json data = nlohmann::json::parse(f);
    return read(data);
  } catch (const std::exception& e) {
    log::fatal("RunReader::read", "Parsing error:", e.what());
  }
  return {};
}

}  // namespace faint
