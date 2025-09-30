#ifndef FAINT_SYST_SYSTEMATICS_H
#define FAINT_SYST_SYSTEMATICS_H

#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"

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

using SystematicTable = std::map<SystematicCategory, std::vector<SystematicDescriptor>>;

inline const std::vector<SystematicDescriptor>& variable_registry_systematics() {
  static const std::vector<SystematicDescriptor> descriptors =
      systematic_list_from_variables();
  return descriptors;
}

inline SystematicTable group_systematics_by_category() {
  SystematicTable grouped;
  for (const auto& descriptor : variable_registry_systematics())
    grouped[descriptor.kind].push_back(descriptor);
  return grouped;
}

inline std::vector<std::string> variable_registry_systematic_names() {
  std::vector<std::string> names;
  names.reserve(variable_registry_systematics().size());
  for (const auto& descriptor : variable_registry_systematics())
    names.push_back(descriptor.name);
  return names;
}

inline ROOT::RDF::RNode nominal_weight_node(ROOT::RDF::RNode df) {
  if (!df.HasColumn("nominal_event_weight"))
    return df.Define("nominal_event_weight", []() { return 1.0; });
  return df;
}

static auto weight_product = [](double nom, double w) {
  const double ww  = (std::isfinite(w) && w > 0.0) ? w : 1.0;
  const double out = nom * ww;
  return (std::isfinite(out) && out >= 0.0) ? out : 1.0;
};

inline ROOT::RDF::RNode universe_weight_node(ROOT::RDF::RNode df,
                                             const SystematicDescriptor& s,
                                             int universe_index) {
  auto node = nominal_weight_node(df);

  switch (s.kind) {
    case SystematicCategory::multisim: {
      if (!node.HasColumn(s.array_column))
        return node.Alias("__w", "nominal_event_weight");
      return node.Define(
          "__w",
          [universe_index](double nom, const auto& v) {
            const double w =
                (universe_index >= 0 &&
                 universe_index < static_cast<int>(v.size()))
                    ? static_cast<double>(v[universe_index])
                    : 1.0;
            return weight_product(nom, w);
          },
          {"nominal_event_weight", s.array_column});
    }
    case SystematicCategory::dual_unisim: {
      const bool use_up = (universe_index == 1); // 0=down,1=up
      const std::string& col = use_up ? s.up_column : s.down_column;
      if (!node.HasColumn(col))
        return node.Alias("__w", "nominal_event_weight");
      return node.Define("__w", weight_product,
                         {"nominal_event_weight", col});
    }
    case SystematicCategory::single_unisim: {
      if (!node.HasColumn(s.single_column))
        return node.Alias("__w", "nominal_event_weight");
      return node.Define("__w", weight_product,
                         {"nominal_event_weight", s.single_column});
    }
  }
  return node;
}

inline ROOT::RDF::RResultPtr<TH1D>
cv_histogram(ROOT::RDF::RNode df,
             const ROOT::RDF::TH1DModel& model,
             const std::string& value_column) {
  auto node = nominal_weight_node(df).Alias("__w", "nominal_event_weight");
  return node.Histo1D(model, value_column, "__w");
}

inline std::vector<ROOT::RDF::RResultPtr<TH1D>>
universe_histograms(ROOT::RDF::RNode df,
                    const SystematicDescriptor& s,
                    const ROOT::RDF::TH1DModel& model,
                    const std::string& value_column) {
  std::vector<ROOT::RDF::RResultPtr<TH1D>> out;
  out.reserve(s.universes);
  for (int u = 0; u < s.universes; ++u) {
    auto node_u = universe_weight_node(df, s, u);
    out.emplace_back(node_u.Histo1D(model, value_column, "__w"));
  }
  return out;
}

inline TMatrixD covariance_matrix_from_histograms(const TH1D* cv_hist,
                                                  const std::vector<const TH1D*>& uni_hists,
                                                  SystematicCategory kind) {
  const int nb = cv_hist->GetNbinsX();
  TMatrixD cov(nb, nb);

  auto binv = [](const TH1D* h, int i) { return h->GetBinContent(i + 1); };

  if (kind == SystematicCategory::multisim) {
    const int n = static_cast<int>(uni_hists.size());
    if (n < 2) return cov;

    std::vector<double> mean(nb, 0.0);
    for (int i = 0; i < nb; ++i) {
      double mi = 0.0;
      for (int u = 0; u < n; ++u) mi += binv(uni_hists[u], i);
      mean[i] = mi / n;
    }
    for (int i = 0; i < nb; ++i) {
      for (int j = 0; j < nb; ++j) {
        double cij = 0.0;
        for (int u = 0; u < n; ++u) {
          const double di = binv(uni_hists[u], i) - mean[i];
          const double dj = binv(uni_hists[u], j) - mean[j];
          cij += di * dj;
        }
        cov[i][j] = cij / (n - 1);
      }
    }
    return cov;
  }

  if (kind == SystematicCategory::single_unisim) {
    if (uni_hists.empty()) return cov;
    const TH1D* h0 = uni_hists.front();
    for (int i = 0; i < nb; ++i) {
      const double di = binv(h0, i) - binv(cv_hist, i);
      for (int j = 0; j < nb; ++j) {
        const double dj = binv(h0, j) - binv(cv_hist, j);
        cov[i][j] = di * dj;
      }
    }
    return cov;
  }

  if (kind == SystematicCategory::dual_unisim) {
    if (uni_hists.size() < 2) return cov;
    const TH1D* hdn = uni_hists[0];
    const TH1D* hup = uni_hists[1];
    for (int i = 0; i < nb; ++i) {
      const double di = 0.5 * (binv(hup, i) - binv(hdn, i)); // (up - down)/2
      for (int j = 0; j < nb; ++j) {
        const double dj = 0.5 * (binv(hup, j) - binv(hdn, j));
        cov[i][j] = di * dj;
      }
    }
    return cov;
  }

  return cov;
}

} // namespace syst
} // namespace faint

#endif // FAINT_SYST_SYSTEMATICS_H
