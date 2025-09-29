#include "faint/RunCatalog.h"

#include <fstream>
#include <stdexcept>
#include <utility>

#include "faint/Log.h"

namespace faint {

namespace {

std::string make_key(const std::string& beam, const std::string& period) {
  return beam + ":" + period;
}

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

  log::fatal("RunCatalog::from_json",
             "Run configuration missing 'beamlines' or 'run_configurations' section");
  throw std::runtime_error("unreachable");
}

}  // namespace

const Run& RunCatalog::get(const std::string& beam,
                           const std::string& period) const {
  const auto key = make_key(beam, period);
  const auto it = runs_.find(key);
  if (it == runs_.end())
    throw std::out_of_range("Run not found: " + key);
  return it->second;
}

RunCatalog RunCatalog::from_json(const nlohmann::json& data) {
  RunCatalog catalog;
  const auto& runs = runs_node(data);
  for (const auto& [beam, periods] : runs.items()) {
    for (const auto& [period, details] : periods.items()) {
      Run run(details, beam, period);
      run.validate();
      const auto key = run.label();
      if (catalog.runs_.count(key) != 0)
        log::fatal("RunCatalog::from_json", "Duplicate run label", key);
      catalog.runs_.emplace(key, std::move(run));
    }
  }
  return catalog;
}

RunCatalog RunCatalog::from_file(const std::string& path) {
  std::ifstream input(path);
  if (!input.is_open())
    log::fatal("RunCatalog::from_file", "Could not open config file", path);
  try {
    const auto data = nlohmann::json::parse(input);
    return from_json(data);
  } catch (const std::exception& e) {
    log::fatal("RunCatalog::from_file", "Parsing error", e.what());
  }
  return {};
}

}  // namespace faint

