#ifndef FAINT_SYSTEMATICS_REGISTRY_H
#define FAINT_SYSTEMATICS_REGISTRY_H

#include <string>
#include <utility>
#include <vector>

#include <faint/Variables.h>

namespace faint {
namespace syst {

// domain noun
enum class SystematicCategory { multisim, single_unisim, dual_unisim };

// domain noun
struct SystematicDescriptor {
  std::string name;         // e.g. "weightsGenie", "RPA"
  SystematicCategory kind;  // multisim / single_unisim / dual_unisim
  int universes = 0;

  // columns (domain nouns)
  std::string array_column; // multisim
  std::string up_column;    // dual_unisim
  std::string down_column;  // dual_unisim
  std::string single_column;// single_unisim
};

// function name: lower_snake_case, noun phrase
inline std::vector<SystematicDescriptor> systematic_list_from_variables() {
  std::vector<SystematicDescriptor> out;

  // dual_unisim knobs
  for (const auto& kv : Variables::knob_var()) {
    SystematicDescriptor s;
    s.name        = kv.first;          // "RPA", "CCMEC", ...
    s.kind        = SystematicCategory::dual_unisim;
    s.universes   = 2;
    s.up_column   = kv.second.first;   // e.g. "knobRPAup"
    s.down_column = kv.second.second;  // e.g. "knobRPAdn"
    out.push_back(std::move(s));
  }

  // multisim arrays
  for (const auto& kv : Variables::multi_uni_var()) {
    SystematicDescriptor s;
    s.name         = kv.first;                     // e.g. "weightsGenie"
    s.kind         = SystematicCategory::multisim;
    s.universes    = static_cast<int>(kv.second);  // e.g. 500
    s.array_column = kv.first;                     // same as name
    out.push_back(std::move(s));
  }

  // single_unisim (optional one-sided)
  {
    SystematicDescriptor s;
    s.name          = Variables::single_knob_var(); // "RootinoFix"
    s.kind          = SystematicCategory::single_unisim;
    s.universes     = 1;
    s.single_column = s.name;
    out.push_back(std::move(s));
  }

  return out;
}

} // namespace syst
} // namespace faint

#endif // FAINT_SYSTEMATICS_REGISTRY_H
