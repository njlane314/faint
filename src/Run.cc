#include "faint/Run.h"

#include <fstream>
#include <set>
#include <stdexcept>
#include <utility>


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

  throw std::runtime_error(
      "RunReader::load_from_json: Run configuration missing 'beamlines' or 'run_configurations' section");
}

}  // namespace

Run::Run(const json& j, std::string bm, std::string pr)
    : beam_mode(std::move(bm)),
      run_period(std::move(pr)),
      nominal_pot(j.value("nominal_pot",
                          j.value("pot_target_wcut_total",
                                  j.value("torb_target_pot_wcut", 0.0)))),
      nominal_triggers(j.value("nominal_triggers",
                               j.value("ext_triggers_total",
                                       j.value("ext_triggers", 0L)))),
      samples(j.at("samples")) {}

std::string Run::key() const { return beam_mode + ":" + run_period; }

std::string Run::label() const { return key(); }

void Run::validate() const {
  if (beam_mode.empty()) {
    throw std::runtime_error("Run::validate: empty beam_mode");
  }
  if (run_period.empty()) {
    throw std::runtime_error("Run::validate: empty run_period");
  }
  if (samples.empty()) {
    throw std::runtime_error("Run::validate: no samples for " + beam_mode +
                             "/" + run_period);
  }

  std::set<std::string> keys;
  for (auto& s : samples) {
    std::string key = s.at("sample_key").get<std::string>();
    if (!keys.insert(key).second) {
      throw std::runtime_error("Run::validate: duplicate sample key: " + key);
    }
  }
}

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

RunReader::RunReader(const std::string& path) {
  std::ifstream f(path);
  if (!f.is_open()) {
    throw std::runtime_error("RunReader::RunReader: Could not open config file: " + path);
  }

  try {
    load_from_json(nlohmann::json::parse(f));
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("RunReader::RunReader: Parsing error: ") +
                             e.what());
  }
}

void RunReader::load_from_json(const nlohmann::json& data) {
  configs_.clear();
  const auto& runs = runs_node(data);
  for (auto const& [beam, periods] : runs.items()) {
    for (auto const& [period, details] : periods.items()) {
      Run rc(details, beam, period);
      rc.validate();
      add(std::move(rc));
    }
  }
}

}  // namespace faint
