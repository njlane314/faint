#pragma once

#include <ROOT/RDataFrame.hxx>
#include <TH1.h>
#include <TH2.h>
#include <THn.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rarexsec {
namespace ps {

struct Axes4 {
    int nb_pmu = 20;
    double min_pmu = 0.1;
    double max_pmu = 5.0;
    int nb_cth_mu = 20;
    double min_cth_mu = -1.0;
    double max_cth_mu = 1.0;
    int nb_log_bg = 20;
    double min_log_bg = -1.0;
    double max_log_bg = 1.3;
    int nb_cth_l = 20;
    double min_cth_l = -1.0;
    double max_cth_l = 1.0;
    const char* histogram_name = "Hgen";
};

class Reweighter {
  public:
    static std::shared_ptr<Reweighter> from_histogram(const THnF& histogram, double laplace = 1.0) {
        auto* clone = static_cast<THnF*>(histogram.Clone((std::string(histogram.GetName()) + "_clone").c_str()));
        clone->SetDirectory(nullptr);
        return std::shared_ptr<Reweighter>(new Reweighter(std::unique_ptr<THnF>(clone), laplace));
    }

    float weight_coords(const std::vector<double>& coords) const {
        if (!histogram_) {
            return 1.0f;
        }
        const int dimensions = histogram_->GetNdimensions();
        if (static_cast<int>(coords.size()) != dimensions) {
            throw std::runtime_error("Reweighter::weight_coords: coords size != histogram ndim");
        }
        std::vector<int> bin_indices(dimensions);
        for (int axis_index = 0; axis_index < dimensions; ++axis_index) {
            auto* axis = histogram_->GetAxis(axis_index);
            int bin = axis->FindBin(coords[axis_index]);
            if (bin < 1) {
                bin = 1;
            }
            if (bin > axis->GetNbins()) {
                bin = axis->GetNbins();
            }
            bin_indices[axis_index] = bin;
        }
        const double bin_content = histogram_->GetBinContent(bin_indices.data());
        const double denominator = bin_content + laplace_;
        return denominator > 0.0 ? static_cast<float>(1.0 / denominator) : 0.0f;
    }

    float weight_four(double muon_momentum,
                      double muon_cos_theta,
                      double lambda_momentum,
                      double lambda_cos_theta) const {
        constexpr double lambda_mass = 1.115683;
        const double log_beta_gamma = std::log10(std::max(1e-12, lambda_momentum / lambda_mass));
        return weight_coords({muon_momentum, muon_cos_theta, log_beta_gamma, lambda_cos_theta});
    }

    int ndim() const noexcept {
        return histogram_ ? histogram_->GetNdimensions() : 0;
    }

    double laplace() const noexcept {
        return laplace_;
    }

  private:
    explicit Reweighter(std::unique_ptr<THnF> histogram, double laplace)
        : histogram_(std::move(histogram)), laplace_(laplace) {}

    std::unique_ptr<THnF> histogram_;
    double laplace_;
};

inline std::unique_ptr<THnF> build_hgen(ROOT::RDF::RNode node,
                                        const Axes4& axes,
                                        std::string_view filter_expr = "is_signal && phi_has_lambda",
                                        std::string_view col_pmu = "phi_pmu",
                                        std::string_view col_cth_mu = "phi_costh_mu",
                                        std::string_view col_p_lambda = "phi_pL",
                                        std::string_view col_cth_lambda = "phi_costh_L") {
    const int dimensions = 4;
    int bins[dimensions] = {axes.nb_pmu, axes.nb_cth_mu, axes.nb_log_bg, axes.nb_cth_l};
    double min_values[dimensions] = {axes.min_pmu, axes.min_cth_mu, axes.min_log_bg, axes.min_cth_l};
    double max_values[dimensions] = {axes.max_pmu, axes.max_cth_mu, axes.max_log_bg, axes.max_cth_l};

    auto master = std::make_unique<THnF>(axes.histogram_name,
                                         "Phase-space; p_{#mu}; cos#theta_{#mu}; log_{10}(#beta#gamma_{#Lambda}); cos#theta_{#Lambda}",
                                         dimensions,
                                         bins,
                                         min_values,
                                         max_values);
    master->Sumw2();

    const unsigned slot_count = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1U;
    std::vector<std::unique_ptr<THnF>> local_histograms(slot_count);

    auto filtered_node = node.Filter(std::string(filter_expr));

    auto action = filtered_node.ForeachSlot(
        [&, bins, min_values, max_values](unsigned slot,
                                          float muon_momentum,
                                          float muon_cos_theta,
                                          float lambda_momentum,
                                          float lambda_cos_theta) {
            auto& histogram = local_histograms[slot];
            if (!histogram) {
                histogram.reset(new THnF((std::string(axes.histogram_name) + "_local_" + std::to_string(slot)).c_str(),
                                         master->GetTitle(),
                                         dimensions,
                                         bins,
                                         min_values,
                                         max_values));
                histogram->Sumw2();
            }
            constexpr double lambda_mass = 1.115683;
            const double log_beta_gamma = std::log10(std::max(1e-12, static_cast<double>(lambda_momentum) / lambda_mass));
            const double coordinates[4] = {static_cast<double>(muon_momentum),
                                           static_cast<double>(muon_cos_theta),
                                           log_beta_gamma,
                                           static_cast<double>(lambda_cos_theta)};
            histogram->Fill(coordinates, 1.0);
        },
        {std::string(col_pmu).c_str(),
         std::string(col_cth_mu).c_str(),
         std::string(col_p_lambda).c_str(),
         std::string(col_cth_lambda).c_str()});

    action.GetValue();

    for (auto& local_histogram : local_histograms) {
        if (local_histogram) {
            master->Add(local_histogram.get());
        }
    }
    return master;
}

inline ROOT::RDF::RNode define_phase_space_weight(ROOT::RDF::RNode node,
                                                  std::shared_ptr<Reweighter> reweighter,
                                                  std::string output_column = "w_ps",
                                                  std::string_view col_has_lambda = "phi_has_lambda",
                                                  std::string_view col_pmu = "phi_pmu",
                                                  std::string_view col_cth_mu = "phi_costh_mu",
                                                  std::string_view col_p_lambda = "phi_pL",
                                                  std::string_view col_cth_lambda = "phi_costh_L") {
    return node.DefineSlot(
        output_column.c_str(),
        [reweighter](unsigned, bool has_lambda, float muon_momentum, float muon_cos_theta, float lambda_momentum, float lambda_cos_theta) -> float {
            if (!has_lambda || muon_momentum <= 0.0f || lambda_momentum <= 0.0f) {
                return 1.0f;
            }
            return reweighter->weight_four(muon_momentum, muon_cos_theta, lambda_momentum, lambda_cos_theta);
        },
        {std::string(col_has_lambda).c_str(),
         std::string(col_pmu).c_str(),
         std::string(col_cth_mu).c_str(),
         std::string(col_p_lambda).c_str(),
         std::string(col_cth_lambda).c_str()});
}

inline ROOT::RDF::RNode attach_weights(ROOT::RDF::RNode node,
                                       const Axes4& axes = Axes4{},
                                       double laplace = 1.0,
                                       std::string_view filter_expr = "is_signal && phi_has_lambda",
                                       std::string output_column = "w_ps",
                                       std::string_view col_has_lambda = "phi_has_lambda",
                                       std::string_view col_pmu = "phi_pmu",
                                       std::string_view col_cth_mu = "phi_costh_mu",
                                       std::string_view col_p_lambda = "phi_pL",
                                       std::string_view col_cth_lambda = "phi_costh_L",
                                       std::shared_ptr<Reweighter>* out_reweighter = nullptr) {
    auto histogram = build_hgen(node, axes, filter_expr, col_pmu, col_cth_mu, col_p_lambda, col_cth_lambda);
    auto reweighter = Reweighter::from_histogram(*histogram, laplace);
    if (out_reweighter) {
        *out_reweighter = reweighter;
    }
    return define_phase_space_weight(node,
                                     std::move(reweighter),
                                     output_column,
                                     col_has_lambda,
                                     col_pmu,
                                     col_cth_mu,
                                     col_p_lambda,
                                     col_cth_lambda);
}

inline double ks_vs_uniform_one_d(const TH1& histogram) {
    std::unique_ptr<TH1> weighted_histogram(static_cast<TH1*>(histogram.Clone("ks_tmp_h")));
    const double integral = weighted_histogram->Integral("width");
    if (integral > 0.0) {
        weighted_histogram->Scale(1.0 / integral);
    }

    std::unique_ptr<TH1> uniform_histogram(static_cast<TH1*>(weighted_histogram->Clone("ks_tmp_uni")));
    uniform_histogram->Reset();
    const int bin_count = uniform_histogram->GetXaxis()->GetNbins();
    const double length = uniform_histogram->GetXaxis()->GetXmax() - uniform_histogram->GetXaxis()->GetXmin();
    const double density = length > 0.0 ? 1.0 / length : 0.0;
    for (int bin = 1; bin <= bin_count; ++bin) {
        uniform_histogram->SetBinContent(bin, density);
    }
    const double uniform_integral = uniform_histogram->Integral("width");
    if (uniform_integral > 0.0) {
        uniform_histogram->Scale(1.0 / uniform_integral);
    }

    return weighted_histogram->KolmogorovTest(uniform_histogram.get());
}

inline std::pair<double, double> chi_squared_constant_two_d(const TH2& histogram) {
    std::unique_ptr<TH2> weighted_histogram(static_cast<TH2*>(histogram.Clone("chi2_tmp")));
    const int bins_x = weighted_histogram->GetXaxis()->GetNbins();
    const int bins_y = weighted_histogram->GetYaxis()->GetNbins();
    const int bin_total = bins_x * bins_y;
    const double total = weighted_histogram->Integral();
    if (total <= 0.0) {
        return {0.0, 1.0};
    }
    const double mean = total / bin_total;

    double chi_squared = 0.0;
    int degrees_of_freedom = 0;
    for (int bin_x = 1; bin_x <= bins_x; ++bin_x) {
        for (int bin_y = 1; bin_y <= bins_y; ++bin_y) {
            const double sum = weighted_histogram->GetBinContent(bin_x, bin_y);
            const double error = weighted_histogram->GetBinError(bin_x, bin_y);
            const double variance = std::max(1e-12, error * error);
            chi_squared += (sum - mean) * (sum - mean) / variance;
            ++degrees_of_freedom;
        }
    }
    degrees_of_freedom = std::max(1, degrees_of_freedom - 1);
    return {chi_squared / degrees_of_freedom, TMath::Prob(chi_squared, degrees_of_freedom)};
}

}
}
