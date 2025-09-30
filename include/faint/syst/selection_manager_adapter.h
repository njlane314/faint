#ifndef FAINT_SYST_SELECTION_MANAGER_ADAPTER_H
#define FAINT_SYST_SELECTION_MANAGER_ADAPTER_H

#include <map>
#include <string>
#include <vector>

#include <faint/syst/systematics_registry.h>

namespace faint {
namespace syst {

using SystematicTable = std::map<SystematicCategory, std::vector<SystematicDescriptor>>;

// function (noun)
inline const std::vector<SystematicDescriptor>& selection_manager_systematics() {
  static const std::vector<SystematicDescriptor> descriptors =
      systematic_list_from_variables();
  return descriptors;
}

// function (noun phrase)
inline SystematicTable group_systematics_by_category() {
  SystematicTable grouped;
  for (const auto& descriptor : selection_manager_systematics())
    grouped[descriptor.kind].push_back(descriptor);
  return grouped;
}

// function (noun phrase)
inline std::vector<std::string> selection_manager_systematic_names() {
  std::vector<std::string> names;
  names.reserve(selection_manager_systematics().size());
  for (const auto& descriptor : selection_manager_systematics())
    names.push_back(descriptor.name);
  return names;
}

} // namespace syst
} // namespace faint

#endif // FAINT_SYST_SELECTION_MANAGER_ADAPTER_H
