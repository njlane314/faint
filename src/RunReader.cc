#include "faint/RunReader.h"

#include <fstream>
#include <stdexcept>
#include <utility>

#include "faint/Log.h"

namespace faint {

namespace {

const nlohmann::json& runs_node(const nlohmann::json& data) {
  const nlohmann::json* node = &data;
  if (node->contains("samples")) {
    node = &node->at("samples");
  }

  if (node->contains("run_configurations")) {
    return node->at("run_configurations");
  }
  if (node->contains("beamlines")) {
    return node->at("beamlines");
  }

  log::fatal("RunReader::from_json",
             "Run configuration missing 'beamlines' or 'run_configurations' section");
  throw std::runtime_error("unreachable");
}

}  // namespace

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

RunReader RunReader::from_json(const nlohmann::json& data) {
  RunReader out;
  const auto& runs = runs_node(data);
  for (auto const& [beam, periods] : runs.items()) {
    for (auto const& [period, details] : periods.items()) {
      Run rc(details, beam, period);
      rc.validate();
      out.add(std::move(rc));
    }
  }
  return out;
}

RunReader RunReader::from_file(const std::string& path) {
  std::ifstream f(path);
  if (!f.is_open())
    log::fatal("RunReader::from_file", "Could not open config file:", path);
  try {
    nlohmann::json data = nlohmann::json::parse(f);
    return from_json(data);
  } catch (const std::exception& e) {
    log::fatal("RunReader::from_file", "Parsing error:", e.what());
  }
  return {};
}

}  // namespace faint
