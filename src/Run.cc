#include "faint/Run.h"

#include <fstream>
#include <set>
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
  if (beam_mode.empty())
    log::fatal("Run::validate", "empty beam_mode");
  if (run_period.empty())
    log::fatal("Run::validate", "empty run_period");
  if (samples.empty())
    log::fatal("Run::validate", "no samples for", beam_mode + "/" + run_period);

  std::set<std::string> keys;
  for (auto& s : samples) {
    std::string key = s.at("sample_key").get<std::string>();
    if (!keys.insert(key).second)
      log::fatal("Run::validate", "duplicate sample key:", key);
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
