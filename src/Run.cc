#include "faint/Run.h"

#include <set>
#include <utility>

#include "faint/utils/Logger.h"

namespace faint {

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

}  // namespace faint
