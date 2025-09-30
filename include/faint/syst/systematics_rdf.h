#ifndef FAINT_SYSTEMATICS_RDF_H
#define FAINT_SYSTEMATICS_RDF_H

#include <cmath>
#include <string>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"

#include <faint/syst/systematics_registry.h>

namespace faint {
namespace syst {

// function (noun phrase)
inline ROOT::RDF::RNode nominal_weight_node(ROOT::RDF::RNode df) {
  if (!df.HasColumn("nominal_event_weight"))
    return df.Define("nominal_event_weight", []() { return 1.0; });
  return df;
}

// helper variable (noun)
static auto weight_product = [](double nom, double w) {
  const double ww  = (std::isfinite(w) && w > 0.0) ? w : 1.0;
  const double out = nom * ww;
  return (std::isfinite(out) && out >= 0.0) ? out : 1.0;
};

// function (noun phrase)
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

// function (noun)
inline ROOT::RDF::RResultPtr<TH1D>
cv_histogram(ROOT::RDF::RNode df,
             const ROOT::RDF::TH1DModel& model,
             const std::string& value_column) {
  auto node = nominal_weight_node(df).Alias("__w", "nominal_event_weight");
  return node.Histo1D(model, value_column, "__w");
}

// function (noun)
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

// function (noun phrase)
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

#endif // FAINT_SYSTEMATICS_RDF_H
