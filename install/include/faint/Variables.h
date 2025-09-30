#ifndef ANALYSIS_VARIABLES_H
#define ANALYSIS_VARIABLES_H

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <faint/Types.h>

namespace faint {

class Variables {
 public:
  using KnobVariations =
      std::unordered_map<std::string, std::pair<std::string, std::string>>;
  using MultiUniverseVars = std::unordered_map<std::string, unsigned>;

  static const KnobVariations& knob_var();

  static const MultiUniverseVars& multi_uni_var();

  static const std::string& single_knob_var();

  static std::vector<std::string> event_var(SampleOrigin origin);

 private:
  static std::unordered_set<std::string> collect_base_vars();
  static void add_mc_vars(std::unordered_set<std::string>& vars);

  static const std::vector<std::string>& base_var();
  static const std::vector<std::string>& truth_var();
  static const std::vector<std::string>& reco_var();
  static const std::vector<std::string>& image_var();
  static const std::vector<std::string>& flash_var();
  static const std::vector<std::string>& energy_var();
  static const std::vector<std::string>& slice_var();
  static const std::vector<std::string>& track_var();
  static const std::vector<std::string>& proc_evt_var();
};

using VariableRegistry = Variables;

}  // namespace faint

#endif
