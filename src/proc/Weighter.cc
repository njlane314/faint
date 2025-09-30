#include "faint/proc/Weighter.h"

#include <cmath>

#include "faint/Log.h"

namespace faint {

using sample::SampleOrigin;

Weighter::Weighter(const nlohmann::json& cfg, double total_run_pot,
                   long total_run_triggers)
    : sample_pot_(cfg.value("pot", 0.0)),
      sample_triggers_(cfg.value("triggers", 0L)),
      total_run_pot_(total_run_pot),
      total_run_triggers_(total_run_triggers) {
  if (sample_pot_ <= 0.0 && sample_triggers_ <= 0L) {
    log::warn("Weighter::Weighter",
              "sample JSON has no or invalid 'pot' or 'triggers'; "
              "base_event_weight will default to 1");
  }
}

ROOT::RDF::RNode Weighter::process(ROOT::RDF::RNode df,
                                   SampleOrigin origin) const {
  ROOT::RDF::RNode node = df;

  if (origin == SampleOrigin::kMonteCarlo || origin == SampleOrigin::kDirt) {
    double scale = 1.0;
    if (sample_pot_ > 0.0 && total_run_pot_ > 0.0)
      scale = total_run_pot_ / sample_pot_;

    node = node.Define("base_event_weight", [scale]() { return scale; });

    node = node.Define(
        "nominal_event_weight",
        [](double w, float w_spline, float w_tune) {
          double out = w;
          if (std::isfinite(w_spline) && w_spline > 0)
            out *= w_spline;
          if (std::isfinite(w_tune) && w_tune > 0)
            out *= w_tune;
          if (!std::isfinite(out) || out < 0)
            return 1.0;
          return out;
        },
        {"base_event_weight", "weightSpline", "weightTune"});

  } else if (origin == SampleOrigin::kExternal) {
    double scale = 1.0;
    if (sample_triggers_ > 0 && total_run_triggers_ > 0) {
      scale = static_cast<double>(total_run_triggers_) /
              static_cast<double>(sample_triggers_);
    }
    node = node.Define("base_event_weight", [scale]() { return scale; });
  }

  if (!node.HasColumn("nominal_event_weight")) {
    if (node.HasColumn("base_event_weight"))
      node = node.Alias("nominal_event_weight", "base_event_weight");
    else
      node = node.Define("nominal_event_weight", []() { return 1.0; });
  }

  return next_ ? next_->process(node, origin) : node;
}

}  // namespace faint
