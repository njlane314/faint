#ifndef FAINT_SELECTION_MANAGER_ADAPTER_H
#define FAINT_SELECTION_MANAGER_ADAPTER_H

#include <string>
#include <vector>

#include "SelectionManager.h"
#include "SystematicsHeader.h"

// function (noun)
inline void systematic_setup(SelectionManager& sel) {
  // Flux
  for (size_t i = 0; i < FluxU_Str.size(); ++i)
    sel.AddSystematic(AllU_SysTypes[sFlux][i], AllU_Universes[sFlux][i], AllU_Str[sFlux][i]);

  // Generator
  for (size_t i = 0; i < GeneratorU_Str.size(); ++i)
    sel.AddSystematic(AllU_SysTypes[sGenerator][i], AllU_Universes[sGenerator][i], AllU_Str[sGenerator][i]);

  // Reinteraction
  for (size_t i = 0; i < ReintU_Str.size(); ++i)
    sel.AddSystematic(AllU_SysTypes[sReint][i], AllU_Universes[sReint][i], AllU_Str[sReint][i]);

  // Misc
  for (size_t i = 0; i < MiscU_Str.size(); ++i)
    sel.AddSystematic(AllU_SysTypes[sMisc][i], AllU_Universes[sMisc][i], AllU_Str[sMisc][i]);

  // Detector
  for (size_t i = 0; i < DetectorU_Str.size(); ++i)
    sel.AddSystematic(AllU_SysTypes[sDetector][i], AllU_Universes[sDetector][i], AllU_Str[sDetector][i]);
}

#endif // FAINT_SELECTION_MANAGER_ADAPTER_H
