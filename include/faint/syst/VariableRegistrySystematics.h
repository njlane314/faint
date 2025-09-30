#ifndef FAINT_SYST_VARIABLE_REGISTRY_SYSTEMATICS_H
#define FAINT_SYST_VARIABLE_REGISTRY_SYSTEMATICS_H

#include <map>
#include <string>
#include <vector>

#include <faint/syst/SystematicsRegistry.h>

namespace faint {
namespace syst {

using SystematicTable = std::map<SystematicCategory, std::vector<SystematicDescriptor>>;

// function (noun)
inline const std::vector<SystematicDescriptor>& variable_registry_systematics() {
  static const std::vector<SystematicDescriptor> descriptors =
      systematic_list_from_variables();
  return descriptors;
}

// function (noun phrase)
inline SystematicTable group_systematics_by_category() {
  SystematicTable grouped;
  for (const auto& descriptor : variable_registry_systematics())
    grouped[descriptor.kind].push_back(descriptor);
  return grouped;
}

// function (noun phrase)
inline std::vector<std::string> variable_registry_systematic_names() {
  std::vector<std::string> names;
  names.reserve(variable_registry_systematics().size());
  for (const auto& descriptor : variable_registry_systematics())
    names.push_back(descriptor.name);
  return names;
}

} // namespace syst
} // namespace faint

#endif // FAINT_SYST_VARIABLE_REGISTRY_SYSTEMATICS_H
