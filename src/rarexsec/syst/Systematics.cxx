#include "rarexsec/syst/Systematics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using MapSD = std::map<std::string, std::vector<double>>;
using rarexsec::plot::TH1DModel;

static std::string expr_column_name(const rarexsec::plot::TH1DModel& spec) {
    const std::string base = !spec.id.empty()
                                 ? spec.id
                                 : (!spec.name.empty() ? spec.name : std::string{"expr"});
    return "_rx_expr_" + rarexsec::plot::Plotter::sanitise(base);
}

ROOT::RDF::RNode with_expr(ROOT::RDF::RNode node, const rarexsec::plot::TH1DModel& spec) {
    if (spec.expr.empty())
        return node;
    const std::string col = expr_column_name(spec);
    return node.Define(col, spec.expr);
}

static std::string expr_var(const rarexsec::plot::TH1DModel& spec) {
    if (spec.expr.empty()) {
        if (!spec.id.empty())
            return spec.id;
        if (!spec.name.empty())
            return spec.name;
        return std::string{"x"};
    }
    return expr_column_name(spec);
}

static std::unique_ptr<TH1D> sum_hists(std::vector<ROOT::RDF::RResultPtr<TH1D>> parts,
                                       const std::string& name) {
    std::unique_ptr<TH1D> total;
    for (auto& rr : parts) {
        const TH1D& h = rr.GetValue();
        if (!total) {
            total.reset(static_cast<TH1D*>(h.Clone(name.c_str())));
            if (total)
                total->SetDirectory(nullptr);
        } else if (total) {
            total->Add(&h);
        }
    }
    return total;
}

std::unique_ptr<TH1D> rarexsec::syst::make_total_mc_hist(const rarexsec::plot::TH1DModel& spec,
                                                         const std::vector<const Entry*>& entries,
                                                         const std::string& suffix) {
    TH1::SetDefaultSumw2(true);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
    parts.reserve(entries.size());
    for (size_t ie = 0; ie < entries.size(); ++ie) {
        const Entry* e = entries[ie];
        if (!e)
            continue;
        auto n0 = selection::apply(e->rnode(), spec.sel, *e);
        auto n1 = with_expr(n0, spec);
        auto var = expr_var(spec);
        parts.push_back(n1.Histo1D(spec.model("_mc_src" + std::to_string(ie) + suffix), var, spec.weight));
    }
    auto hist = sum_hists(std::move(parts), spec.id + suffix);
    if (!hist) {
        const std::string name = rarexsec::plot::Plotter::sanitise(spec.id + suffix);
        const std::string title = spec.title.empty() ? spec.id : spec.title;
        auto empty = std::make_unique<TH1D>(name.c_str(), title.c_str(), spec.nbins, spec.xmin, spec.xmax);
        empty->SetDirectory(nullptr);
        return empty;
    }
    return hist;
}

std::unique_ptr<TH1D> rarexsec::syst::make_total_mc_hist_detvar(const rarexsec::plot::TH1DModel& spec,
                                                                const std::vector<const Entry*>& entries,
                                                                const std::string& tag,
                                                                const std::string& suffix) {
    TH1::SetDefaultSumw2(true);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
    parts.reserve(entries.size());
    for (size_t ie = 0; ie < entries.size(); ++ie) {
        const Entry* e = entries[ie];
        if (!e)
            continue;
        const Frame* dv = e->detvar(tag);
        if (!dv)
            continue;
        auto n0 = selection::apply(dv->rnode(), spec.sel, *e);
        auto n1 = with_expr(n0, spec);
        auto var = expr_var(spec);
        parts.push_back(n1.Histo1D(spec.model("_mc_detvar_" + tag + "_src" + std::to_string(ie) + suffix),
                                   var, spec.weight));
    }
    auto hist = sum_hists(std::move(parts), spec.id + suffix);
    if (!hist) {
        const std::string name = rarexsec::plot::Plotter::sanitise(spec.id + suffix);
        const std::string title = spec.title.empty() ? spec.id : spec.title;
        auto empty = std::make_unique<TH1D>(name.c_str(), title.c_str(), spec.nbins, spec.xmin, spec.xmax);
        empty->SetDirectory(nullptr);
        return empty;
    }
    return hist;
}

TMatrixDSym rarexsec::syst::mc_stat_covariance(const TH1D& hist) {
    const int nb = hist.GetNbinsX();
    TMatrixDSym C(nb);
    for (int i = 1; i <= nb; ++i) {
        const double e = hist.GetBinError(i);
        const double v = (std::isfinite(e) && e > 0.0) ? e * e : 0.0;
        C(i - 1, i - 1) = v;
    }
    return C;
}

TMatrixDSym rarexsec::syst::sample_covariance(const TH1D& nominal,
                                              const std::vector<std::unique_ptr<TH1D>>& universes) {
    const int nb = nominal.GetNbinsX();
    std::vector<std::vector<double>> deltas;
    deltas.reserve(universes.size());
    for (const auto& uptr : universes) {
        if (!uptr)
            continue;
        std::vector<double> diff(nb, 0.0);
        for (int i = 0; i < nb; ++i) {
            diff[i] = uptr->GetBinContent(i + 1) - nominal.GetBinContent(i + 1);
        }
        deltas.emplace_back(std::move(diff));
    }
    const int N = static_cast<int>(deltas.size());
    TMatrixDSym C(nb);
    if (N <= 1)
        return C;
    for (int i = 0; i < nb; ++i) {
        for (int j = i; j < nb; ++j) {
            long double s = 0.0L;
            for (const auto& diff : deltas)
                s += static_cast<long double>(diff[i]) * diff[j];
            const double cij = static_cast<double>(s / (N - 1));
            C(i, j) = C(j, i) = cij;
        }
    }
    return C;
}

TMatrixDSym rarexsec::syst::hessian_covariance(const TH1D& nominal,
                                               const TH1D& up,
                                               const TH1D& down) {
    const int nb = nominal.GetNbinsX();
    TMatrixDSym C(nb);
    for (int i = 0; i < nb; ++i) {
        const double dpi = up.GetBinContent(i + 1) - nominal.GetBinContent(i + 1);
        const double dmi = down.GetBinContent(i + 1) - nominal.GetBinContent(i + 1);
        for (int j = i; j < nb; ++j) {
            const double dpj = up.GetBinContent(j + 1) - nominal.GetBinContent(j + 1);
            const double dmj = down.GetBinContent(j + 1) - nominal.GetBinContent(j + 1);
            const double cij = 0.5 * (dpi * dpj + dmi * dmj);
            C(i, j) = C(j, i) = cij;
        }
    }
    return C;
}

TMatrixDSym rarexsec::syst::sum(const std::vector<const TMatrixDSym*>& pieces) {
    if (pieces.empty())
        return TMatrixDSym(0);
    const TMatrixDSym* first = nullptr;
    for (const auto* p : pieces) {
        if (p) {
            first = p;
            break;
        }
    }
    if (!first)
        return TMatrixDSym(0);
    TMatrixDSym C(first->GetNrows());
    for (const auto* p : pieces) {
        if (!p)
            continue;
        C += *p;
    }
    return C;
}

std::unique_ptr<TH1D> rarexsec::syst::make_total_mc_hist_weight_universe_ushort(
    const TH1DModel& spec, const std::vector<const Entry*>& mc,
    const std::string& weights_branch, int k, const std::string& suffix,
    const std::string& cv_branch, double us_scale) {

    TH1::SetDefaultSumw2(true);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
    parts.reserve(mc.size());

    for (size_t ie = 0; ie < mc.size(); ++ie) {
        const Entry* e = mc[ie];
        if (!e)
            continue;

        auto n0 = selection::apply(e->rnode(), spec.sel, *e);
        auto n1 = with_expr(n0, spec);
        auto var = expr_var(spec);

        const std::string col = "_w_us_univ_" + std::to_string(k) + "_src" + std::to_string(ie);
        if (cv_branch.empty()) {
            auto n2 = n1.Define(
                col,
                [k, us_scale](const ROOT::RVec<unsigned short>& v, double w_nom) {
                    double wk = 1.0;
                    if (k >= 0 && k < (int)v.size())
                        wk = static_cast<double>(v[k]) * us_scale;
                    const double out = w_nom * wk;
                    return std::isfinite(out) && out > 0.0 ? out : 0.0;
                },
                {weights_branch, spec.weight});
            parts.push_back(n2.Histo1D(spec.model("_mc_univ_us_" + std::to_string(k) + "_src" + std::to_string(ie) + suffix),
                                       var, col));
        } else {
            auto n2 = n1.Define(
                col,
                [k, us_scale](const ROOT::RVec<unsigned short>& v, double w_nom, double w_cv) {
                    double wk = 1.0;
                    if (k >= 0 && k < (int)v.size())
                        wk = static_cast<double>(v[k]) * us_scale;
                    const double out = w_nom * w_cv * wk;
                    return std::isfinite(out) && out > 0.0 ? out : 0.0;
                },
                {weights_branch, spec.weight, cv_branch});
            parts.push_back(n2.Histo1D(spec.model("_mc_univ_us_" + std::to_string(k) + "_src" + std::to_string(ie) + suffix),
                                       var, col));
        }
    }
    return sum_hists(std::move(parts), spec.id + suffix);
}

TMatrixDSym rarexsec::syst::cov_from_weight_vector_ushort(
    const TH1DModel& spec, const std::vector<const Entry*>& mc,
    const std::string& weights_branch, int nuniv,
    const std::string& cv_branch, double us_scale) {

    if (nuniv <= 0)
        return TMatrixDSym(0);
    auto H0 = rarexsec::syst::make_total_mc_hist(spec, mc, "_nom");
    std::vector<std::unique_ptr<TH1D>> universes;
    universes.reserve(nuniv);
    for (int k = 0; k < nuniv; ++k) {
        universes.emplace_back(
            rarexsec::syst::make_total_mc_hist_weight_universe_ushort(spec, mc, weights_branch, k,
                                                                      "_us_" + std::to_string(k),
                                                                      cv_branch, us_scale));
    }
    return rarexsec::syst::sample_covariance(*H0, universes);
}

std::unique_ptr<TH1D> rarexsec::syst::make_total_mc_hist_weight_universe_map(
    const TH1DModel& spec, const std::vector<const Entry*>& mc,
    const std::string& map_branch, const std::string& key, int k,
    const std::string& suffix, const std::string& cv_branch) {

    TH1::SetDefaultSumw2(true);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
    parts.reserve(mc.size());

    for (size_t ie = 0; ie < mc.size(); ++ie) {
        const Entry* e = mc[ie];
        if (!e)
            continue;

        auto n0 = selection::apply(e->rnode(), spec.sel, *e);
        auto n1 = with_expr(n0, spec);
        auto var = expr_var(spec);

        const std::string col = "_w_map_univ_" + std::to_string(k) + "_src" + std::to_string(ie);
        if (cv_branch.empty()) {
            auto n2 = n1.Define(
                col,
                [k, key](const MapSD& m, double w_nom) {
                    double wk = 1.0;
                    auto it = m.find(key);
                    if (it != m.end() && k >= 0 && k < (int)it->second.size())
                        wk = it->second[k];
                    const double out = w_nom * wk;
                    return std::isfinite(out) && out > 0.0 ? out : 0.0;
                },
                {map_branch, spec.weight});
            parts.push_back(n2.Histo1D(spec.model("_mc_univ_map_" + std::to_string(k) + "_src" + std::to_string(ie) + suffix),
                                       var, col));
        } else {
            auto n2 = n1.Define(
                col,
                [k, key](const MapSD& m, double w_nom, double w_cv) {
                    double wk = 1.0;
                    auto it = m.find(key);
                    if (it != m.end() && k >= 0 && k < (int)it->second.size())
                        wk = it->second[k];
                    const double out = w_nom * w_cv * wk;
                    return std::isfinite(out) && out > 0.0 ? out : 0.0;
                },
                {map_branch, spec.weight, cv_branch});
            parts.push_back(n2.Histo1D(spec.model("_mc_univ_map_" + std::to_string(k) + "_src" + std::to_string(ie) + suffix),
                                       var, col));
        }
    }
    return sum_hists(std::move(parts), spec.id + suffix);
}

TMatrixDSym rarexsec::syst::cov_from_map_weight_vector(
    const TH1DModel& spec, const std::vector<const Entry*>& mc,
    const std::string& map_branch, const std::string& key, int nuniv,
    const std::string& cv_branch) {

    if (nuniv <= 0)
        return TMatrixDSym(0);
    auto H0 = rarexsec::syst::make_total_mc_hist(spec, mc, "_nom");
    std::vector<std::unique_ptr<TH1D>> universes;
    universes.reserve(nuniv);
    for (int k = 0; k < nuniv; ++k) {
        universes.emplace_back(
            rarexsec::syst::make_total_mc_hist_weight_universe_map(spec, mc, map_branch, key, k,
                                                                   "_map_" + std::to_string(k), cv_branch));
    }
    return rarexsec::syst::sample_covariance(*H0, universes);
}

TMatrixDSym rarexsec::syst::cov_from_detvar_pairs(
    const TH1DModel& spec,
    const std::vector<const Entry*>& mc,
    const std::vector<std::pair<std::string, std::string>>& tag_pairs) {

    if (tag_pairs.empty())
        return TMatrixDSym(0);
    auto H0 = rarexsec::syst::make_total_mc_hist(spec, mc, "_nom");
    if (!H0)
        throw std::runtime_error("cov_from_detvar_pairs: failed to build nominal histogram");

    const int nb = H0->GetNbinsX();
    TMatrixDSym C(nb);

    for (const auto& pr : tag_pairs) {
        const auto& up = pr.first;
        const auto& down = pr.second;
        auto Hup = rarexsec::syst::make_total_mc_hist_detvar(spec, mc, up, "_up");
        auto Hdown = rarexsec::syst::make_total_mc_hist_detvar(spec, mc, down, "_down");
        if (!Hup || !Hdown) {
            throw std::runtime_error(
                "cov_from_detvar_pairs: missing detvar hist for tags '" + up + "', '" + down + "'");
        }
        C += rarexsec::syst::hessian_covariance(*H0, *Hup, *Hdown);
    }
    return C;
}

TMatrixDSym rarexsec::syst::cov_from_detvar_unisims(
    const TH1DModel& spec,
    const std::vector<const Entry*>& mc,
    const std::vector<std::string>& tags) {

    if (tags.empty())
        return TMatrixDSym(0);
    auto H0 = rarexsec::syst::make_total_mc_hist(spec, mc, "_nom");
    if (!H0)
        throw std::runtime_error("cov_from_detvar_unisims: failed to build nominal histogram");

    std::vector<std::unique_ptr<TH1D>> universes;
    universes.reserve(tags.size());
    for (const auto& t : tags) {
        auto Ht = rarexsec::syst::make_total_mc_hist_detvar(spec, mc, t, "_var");
        if (!Ht) {
            throw std::runtime_error(
                "cov_from_detvar_unisims: missing detvar hist for tag '" + t + "'");
        }
        universes.emplace_back(std::move(Ht));
    }
    return rarexsec::syst::sample_covariance(*H0, universes);
}

TMatrixDSym rarexsec::syst::block_cov_from_weight_vector_ushort_scaled(
    const TH1DModel& specA, const std::vector<const Entry*>& A,
    const TH1DModel& specB, const std::vector<const Entry*>& B,
    const std::string& weights_branch, int nuniv,
    const std::string& cv_branch, double us_scale) {

    auto H0A = rarexsec::syst::make_total_mc_hist(specA, A, "_A_nom");
    auto H0B = rarexsec::syst::make_total_mc_hist(specB, B, "_B_nom");
    const int nA = H0A ? H0A->GetNbinsX() : 0;
    const int nB = H0B ? H0B->GetNbinsX() : 0;
    TMatrixDSym C(nA + nB);
    if (!H0A || !H0B || nuniv <= 0)
        return C;

    const int ddof = RAREXSEC_MULTISIM_DDOF;
    for (int k = 0; k < nuniv; ++k) {
        auto Au = rarexsec::syst::make_total_mc_hist_weight_universe_ushort(specA, A, weights_branch, k, "_A", cv_branch, us_scale);
        auto Bu = rarexsec::syst::make_total_mc_hist_weight_universe_ushort(specB, B, weights_branch, k, "_B", cv_branch, us_scale);

        std::vector<double> d(nA + nB);
        for (int i = 0; i < nA; ++i)
            d[i] = Au->GetBinContent(i + 1) - H0A->GetBinContent(i + 1);
        for (int j = 0; j < nB; ++j)
            d[nA + j] = Bu->GetBinContent(j + 1) - H0B->GetBinContent(j + 1);

        for (int p = 0; p < nA + nB; ++p)
            for (int q = p; q < nA + nB; ++q)
                C(p, q) += d[p] * d[q];
    }
    C *= 1.0 / std::max(1, nuniv - ddof);
    for (int i = 0; i < nA + nB; ++i)
        for (int j = i + 1; j < nA + nB; ++j)
            C(i, j) = C(j, i) = 0.5 * (C(i, j) + C(j, i));
    return C;
}

TMatrixDSym rarexsec::syst::block_cov_from_map_weight_vector(
    const TH1DModel& specA, const std::vector<const Entry*>& A,
    const TH1DModel& specB, const std::vector<const Entry*>& B,
    const std::string& map_branch, const std::string& key, int nuniv,
    const std::string& cv_branch) {

    auto H0A = rarexsec::syst::make_total_mc_hist(specA, A, "_A_nom");
    auto H0B = rarexsec::syst::make_total_mc_hist(specB, B, "_B_nom");
    const int nA = H0A ? H0A->GetNbinsX() : 0;
    const int nB = H0B ? H0B->GetNbinsX() : 0;
    TMatrixDSym C(nA + nB);
    if (!H0A || !H0B || nuniv <= 0)
        return C;

    const int ddof = RAREXSEC_MULTISIM_DDOF;
    for (int k = 0; k < nuniv; ++k) {
        auto Au = rarexsec::syst::make_total_mc_hist_weight_universe_map(specA, A, map_branch, key, k, "_A", cv_branch);
        auto Bu = rarexsec::syst::make_total_mc_hist_weight_universe_map(specB, B, map_branch, key, k, "_B", cv_branch);

        std::vector<double> d(nA + nB);
        for (int i = 0; i < nA; ++i)
            d[i] = Au->GetBinContent(i + 1) - H0A->GetBinContent(i + 1);
        for (int j = 0; j < nB; ++j)
            d[nA + j] = Bu->GetBinContent(j + 1) - H0B->GetBinContent(j + 1);

        for (int p = 0; p < nA + nB; ++p)
            for (int q = p; q < nA + nB; ++q)
                C(p, q) += d[p] * d[q];
    }
    C *= 1.0 / std::max(1, nuniv - ddof);
    for (int i = 0; i < nA + nB; ++i)
        for (int j = i + 1; j < nA + nB; ++j)
            C(i, j) = C(j, i) = 0.5 * (C(i, j) + C(j, i));
    return C;
}

TMatrixDSym rarexsec::syst::block_cov_from_ud_ushort(
    const TH1DModel& specA, const std::vector<const Entry*>& A,
    const TH1DModel& specB, const std::vector<const Entry*>& B,
    const std::string& up_branch, const std::string& dn_branch, int knob_index,
    double us_scale, const std::string& cv_branch) {

    auto H0A = rarexsec::syst::make_total_mc_hist(specA, A, "_A_nom");
    auto H0B = rarexsec::syst::make_total_mc_hist(specB, B, "_B_nom");
    if (!H0A || !H0B)
        return TMatrixDSym(0);

    auto apply_ud = [&](const TH1DModel& spec, const std::vector<const Entry*>& mc,
                        const std::string& branch, const char* tag) {
        std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
        for (size_t ie = 0; ie < mc.size(); ++ie) {
            const Entry* e = mc[ie];
            if (!e)
                continue;
            auto n0 = selection::apply(e->rnode(), spec.sel, *e);
            auto n1 = with_expr(n0, spec);
            auto var = expr_var(spec);
            const std::string col = std::string("_w_ud_") + tag + "_" + std::to_string(knob_index) + "_src" + std::to_string(ie);

            if (cv_branch.empty()) {
                auto n2 = n1.Define(
                    col,
                    [knob_index, us_scale](const ROOT::RVec<unsigned short>& v, double w_nom) {
                        double wk = 1.0;
                        if (knob_index >= 0 && knob_index < (int)v.size())
                            wk = static_cast<double>(v[knob_index]) * us_scale;
                        const double out = w_nom * wk;
                        return std::isfinite(out) && out > 0.0 ? out : 0.0;
                    },
                    {branch, spec.weight});
                parts.push_back(n2.Histo1D(spec.model(std::string("_mc_ud_") + tag + "_src" + std::to_string(ie)), var, col));
            } else {
                auto n2 = n1.Define(
                    col,
                    [knob_index, us_scale](const ROOT::RVec<unsigned short>& v, double w_nom, double w_cv) {
                        double wk = 1.0;
                        if (knob_index >= 0 && knob_index < (int)v.size())
                            wk = static_cast<double>(v[knob_index]) * us_scale;
                        const double out = w_nom * w_cv * wk;
                        return std::isfinite(out) && out > 0.0 ? out : 0.0;
                    },
                    {branch, spec.weight, cv_branch});
                parts.push_back(n2.Histo1D(spec.model(std::string("_mc_ud_") + tag + "_src" + std::to_string(ie)), var, col));
            }
        }
        return sum_hists(std::move(parts), spec.id + "_" + tag);
    };

    auto HupA = apply_ud(specA, A, up_branch, "upA");
    auto HdnA = apply_ud(specA, A, dn_branch, "dnA");
    auto HupB = apply_ud(specB, B, up_branch, "upB");
    auto HdnB = apply_ud(specB, B, dn_branch, "dnB");

    const int nA = H0A->GetNbinsX();
    const int nB = H0B->GetNbinsX();
    TMatrixDSym C(nA + nB);
    auto hess_cat = [&](const TH1D& H0, const TH1D& Hup, const TH1D& Hdn, int off) {
        const int nb = H0.GetNbinsX();
        for (int i = 1; i <= nb; ++i) {
            const double dpi = Hup.GetBinContent(i) - H0.GetBinContent(i);
            const double dmi = Hdn.GetBinContent(i) - H0.GetBinContent(i);
            for (int j = i; j <= nb; ++j) {
                const double dpj = Hup.GetBinContent(j) - H0.GetBinContent(j);
                const double dmj = Hdn.GetBinContent(j) - H0.GetBinContent(j);
                const double cij = 0.5 * (dpi * dpj + dmi * dmj);
                C(off + i - 1, off + j - 1) += cij;
                C(off + j - 1, off + i - 1) += cij;
            }
        }
    };
    hess_cat(*H0A, *HupA, *HdnA, 0);
    hess_cat(*H0B, *HupB, *HdnB, nA);
    return C;
}

TMatrixDSym rarexsec::syst::block_cov_from_detvar_pairs(
    const TH1DModel& specA, const std::vector<const Entry*>& A,
    const TH1DModel& specB, const std::vector<const Entry*>& B,
    const std::vector<std::pair<std::string, std::string>>& tag_pairs) {

    if (tag_pairs.empty())
        return TMatrixDSym(0);

    auto H0A = rarexsec::syst::make_total_mc_hist(specA, A, "_A_nom");
    auto H0B = rarexsec::syst::make_total_mc_hist(specB, B, "_B_nom");
    if (!H0A || !H0B)
        throw std::runtime_error("block_cov_from_detvar_pairs: failed to build nominal hist(s)");

    const int nA = H0A->GetNbinsX();
    const int nB = H0B->GetNbinsX();
    TMatrixDSym C(nA + nB);

    for (const auto& pr : tag_pairs) {
        const auto& up = pr.first;
        const auto& down = pr.second;

        auto HupA = rarexsec::syst::make_total_mc_hist_detvar(specA, A, up, "_A_up");
        auto HdnA = rarexsec::syst::make_total_mc_hist_detvar(specA, A, down, "_A_dn");
        auto HupB = rarexsec::syst::make_total_mc_hist_detvar(specB, B, up, "_B_up");
        auto HdnB = rarexsec::syst::make_total_mc_hist_detvar(specB, B, down, "_B_dn");
        if (!HupA || !HdnA || !HupB || !HdnB) {
            throw std::runtime_error(
                "block_cov_from_detvar_pairs: missing detvar hist(s) for tags '" + up + "', '" + down + "'");
        }

        std::vector<double> dplus(nA + nB), dminus(nA + nB);
        for (int i = 0; i < nA; ++i) {
            dplus[i] = HupA->GetBinContent(i + 1) - H0A->GetBinContent(i + 1);
            dminus[i] = HdnA->GetBinContent(i + 1) - H0A->GetBinContent(i + 1);
        }
        for (int j = 0; j < nB; ++j) {
            dplus[nA + j] = HupB->GetBinContent(j + 1) - H0B->GetBinContent(j + 1);
            dminus[nA + j] = HdnB->GetBinContent(j + 1) - H0B->GetBinContent(j + 1);
        }

        for (int p = 0; p < nA + nB; ++p) {
            for (int q = p; q < nA + nB; ++q) {
                const double add = 0.5 * (dplus[p] * dplus[q] + dminus[p] * dminus[q]);
                C(p, q) = C(q, p) = C(p, q) + add;
            }
        }
    }
    return C;
}

TMatrixDSym rarexsec::syst::block_diag_stat(const TH1D& A, const TH1D& B) {
    const int nA = A.GetNbinsX();
    const int nB = B.GetNbinsX();
    TMatrixDSym C(nA + nB);
    for (int i = 1; i <= nA; ++i) {
        const double e = A.GetBinError(i);
        C(i - 1, i - 1) = (std::isfinite(e) && e > 0.0) ? e * e : 0.0;
    }
    for (int j = 1; j <= nB; ++j) {
        const double e = B.GetBinError(j);
        C(nA + j - 1, nA + j - 1) = (std::isfinite(e) && e > 0.0) ? e * e : 0.0;
    }
    return C;
}

TMatrixDSym rarexsec::syst::pot_cov_block(const TH1D& A, const TH1D& B, double frac_pot) {
    const int nA = A.GetNbinsX();
    const int nB = B.GetNbinsX();
    const int n = nA + nB;
    std::vector<double> v(n);
    for (int i = 1; i <= nA; ++i)
        v[i - 1] = A.GetBinContent(i);
    for (int j = 1; j <= nB; ++j)
        v[nA + j - 1] = B.GetBinContent(j);
    TMatrixDSym C(n);
    const double s2 = frac_pot * frac_pot;
    for (int p = 0; p < n; ++p)
        for (int q = p; q < n; ++q)
            C(p, q) = C(q, p) = s2 * v[p] * v[q];
    return C;
}

std::unique_ptr<TH1D> rarexsec::syst::sum_same_binning(const TH1D& A, const TH1D& B, const std::string& name) {
    if (A.GetNbinsX() != B.GetNbinsX())
        throw std::runtime_error("sum_same_binning: bin mismatch");
    if (A.GetXaxis()->GetXmin() != B.GetXaxis()->GetXmin() ||
        A.GetXaxis()->GetXmax() != B.GetXaxis()->GetXmax())
        throw std::runtime_error("sum_same_binning: axis range mismatch");

    std::unique_ptr<TH1D> H(static_cast<TH1D*>(A.Clone(name.c_str())));
    H->SetDirectory(nullptr);
    for (int i = 1; i <= A.GetNbinsX(); ++i) {
        const double y = A.GetBinContent(i) + B.GetBinContent(i);
        H->SetBinContent(i, y);
        H->SetBinError(i, 0.0);
    }
    return H;
}

TMatrixDSym rarexsec::syst::sum_covariance_block_same_binning(const TMatrixDSym& C_block, int nA, int nB) {
    if (nA <= 0 || nB <= 0 || C_block.GetNrows() != nA + nB)
        throw std::runtime_error("sum_covariance_block_same_binning: size mismatch");
    const int n = nA;
    if (nB != nA)
        throw std::runtime_error("sum_covariance_block_same_binning: A and B must have the same binning");

    TMatrixDSym Csum(n);
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double cij =
                C_block(i, j) +
                C_block(i, nA + j) +
                C_block(nA + i, j) +
                C_block(nA + i, nA + j);
            Csum(i, j) = Csum(j, i) = cij;
        }
    }
    return Csum;
}
