#include "rarexsec/syst/Systematics.hh"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <ROOT/RDFHelpers.hxx>
#include <TH1.h>

using rarexsec::plot::H1Spec;
using rarexsec::selection::apply;

namespace {

std::string expr_var(const H1Spec& spec) {
    return spec.expr.empty() ? spec.id : "_rx_expr_";
}

ROOT::RDF::RNode with_expr(ROOT::RDF::RNode n, const H1Spec& spec) {
    return spec.expr.empty() ? n : n.Define("_rx_expr_", spec.expr);
}

std::unique_ptr<TH1D> sum_hists(const std::vector<ROOT::RDF::RResultPtr<TH1D>>& parts,
                                const std::string& name) {
    std::unique_ptr<TH1D> total;
    for (auto rr : parts) {
        const TH1D& h = rr.GetValue();
        if (!total) {
            total.reset(static_cast<TH1D*>(h.Clone(name.c_str())));
            total->SetDirectory(nullptr);
        } else {
            total->Add(&h);
        }
    }
    return total;
}

}

namespace rarexsec::syst {

TMatrixDSym mc_stat_covariance(const TH1D& h) {
    const int nb = h.GetNbinsX();
    TMatrixDSym C(nb);
    for (int i = 1; i <= nb; ++i) {
        const double e2 = h.GetBinError(i) * h.GetBinError(i);
        C(i-1,i-1) = std::max(0.0, e2);
    }
    return C;
}

TMatrixDSym sample_covariance(const TH1D& nominal,
                              const std::vector<std::unique_ptr<TH1D>>& universes) {
    const int nb = nominal.GetNbinsX();
    TMatrixDSym C(nb);
    const int N = static_cast<int>(universes.size());
    if (N <= 1) return C;
    std::vector<std::vector<double>> deltas(N, std::vector<double>(nb, 0.0));
    for (int u = 0; u < N; ++u) {
        const TH1D& Hu = *universes[u];
        for (int i = 1; i <= nb; ++i)
            deltas[u][i-1] = Hu.GetBinContent(i) - nominal.GetBinContent(i);
    }
    for (int i = 0; i < nb; ++i) {
        for (int j = i; j < nb; ++j) {
            long double s = 0.0L;
            for (int u = 0; u < N; ++u) s += deltas[u][i] * deltas[u][j];
            const double cij = (N > 1) ? static_cast<double>(s / (N - 1)) : 0.0;
            C(i,j) = C(j,i) = cij;
        }
    }
    return C;
}

TMatrixDSym hessian_covariance(const TH1D& nominal,
                               const TH1D& plus,
                               const TH1D& minus) {
    const int nb = nominal.GetNbinsX();
    TMatrixDSym C(nb);
    for (int i = 1; i <= nb; ++i) {
        const double dplus  = plus.GetBinContent(i)  - nominal.GetBinContent(i);
        const double dminus = minus.GetBinContent(i) - nominal.GetBinContent(i);
        for (int j = 1; j <= nb; ++j) {
            const double eplus  = plus.GetBinContent(j)  - nominal.GetBinContent(j);
            const double eminus = minus.GetBinContent(j) - nominal.GetBinContent(j);
            C(i-1,j-1) = 0.5 * (dplus*eplus + dminus*eminus);
        }
    }
    return C;
}

TMatrixDSym sum(const std::vector<const TMatrixDSym*>& terms) {
    if (terms.empty()) return TMatrixDSym(0);
    const int nb = terms.front()->GetNrows();
    TMatrixDSym C(nb);
    for (auto* t : terms) {
        if (!t) continue;
        if (t->GetNrows() != nb) throw std::runtime_error("Covariance size mismatch in sum()");
        C += *t;
    }
    return C;
}

TMatrixDSym shape_only(const TMatrixDSym& cov, const TH1D& nominal) {
    const int nb = cov.GetNrows();
    TMatrixDSym C = cov;
    std::vector<double> v(nb, 0.0);
    double norm = 0.0;
    for (int i = 1; i <= nb; ++i) { v[i-1] = nominal.GetBinContent(i); norm += v[i-1]*v[i-1]; }
    if (norm <= 0.0) return C;
    for (double& x : v) x /= std::sqrt(norm);
    std::vector<double> u(nb, 0.0);
    long double alpha = 0.0L;
    for (int i = 0; i < nb; ++i) {
        long double ui = 0.0L;
        for (int j = 0; j < nb; ++j) ui += C(i,j) * v[j];
        u[i] = static_cast<double>(ui);
        alpha += v[i] * u[i];
    }
    if (alpha <= 0.0) return C;
    for (int i = 0; i < nb; ++i)
        for (int j = i; j < nb; ++j)
            C(i,j) = C(j,i) = C(i,j) - (u[i]*u[j]/static_cast<double>(alpha));
    return C;
}

std::unique_ptr<TH1D> make_total_mc_hist(const H1Spec& spec,
                                         const std::vector<const Entry*>& mc,
                                         const std::string& suffix) {
    TH1::SetDefaultSumw2(true);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
    for (size_t ie = 0; ie < mc.size(); ++ie) {
        const Entry* e = mc[ie];
        if (!e) continue;
        auto n0 = apply(e->rnode(), spec.sel, *e);
        auto n  = with_expr(n0, spec);
        auto var = expr_var(spec);
        parts.push_back(n.Histo1D(spec.model("_mc_sum_src"+std::to_string(ie)+suffix), var, spec.weight));
    }
    return sum_hists(parts, spec.id + suffix);
}

std::unique_ptr<TH1D> make_total_mc_hist_weight_universe(const H1Spec& spec,
                                                         const std::vector<const Entry*>& mc,
                                                         const std::string& weights_branch,
                                                         int k,
                                                         const std::string& suffix) {
    TH1::SetDefaultSumw2(true);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
    for (size_t ie = 0; ie < mc.size(); ++ie) {
        const Entry* e = mc[ie];
        if (!e) continue;
        auto n0 = apply(e->rnode(), spec.sel, *e);
        auto n1 = with_expr(n0, spec);
        auto var = expr_var(spec);
        const std::string col = "_w_univ_" + std::to_string(k) + "_src" + std::to_string(ie);
        auto n2 = n1.Define(col, [k](const ROOT::RVec<float>& v, float w_nom) {
                                float wk = (k < (int)v.size() && std::isfinite(v[k]) && v[k] > 0.f) ? v[k] : 1.f;
                                float out = w_nom * wk;
                                if (!std::isfinite(out) || out < 0.f) out = 0.f;
                                return out;
                             }, {weights_branch, spec.weight});
        parts.push_back(n2.Histo1D(spec.model("_mc_univ_"+std::to_string(k)+"_src"+std::to_string(ie)+suffix),
                                   var, col));
    }
    return sum_hists(parts, spec.id + suffix);
}

std::unique_ptr<TH1D> make_total_mc_hist_detvar(const H1Spec& spec,
                                                const std::vector<const Entry*>& mc,
                                                const std::string& tag,
                                                const std::string& suffix) {
    TH1::SetDefaultSumw2(true);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> parts;
    for (size_t ie = 0; ie < mc.size(); ++ie) {
        const Entry* e = mc[ie];
        if (!e) continue;
        const auto* dv = e->detvar(tag);
        if (!dv || !dv->node) continue;
        auto n0 = apply(dv->rnode(), spec.sel, *e);
        auto n1 = with_expr(n0, spec);
        auto var = expr_var(spec);
        parts.push_back(n1.Histo1D(spec.model("_mc_detvar_"+tag+"_src"+std::to_string(ie)+suffix),
                                   var, spec.weight));
    }
    return sum_hists(parts, spec.id + suffix);
}

TMatrixDSym cov_from_weight_vector(const H1Spec& spec,
                                   const std::vector<const Entry*>& mc,
                                   const std::string& weights_branch,
                                   int nuniv) {
    if (nuniv <= 0) return TMatrixDSym(0);
    auto H0 = make_total_mc_hist(spec, mc, "_nom");
    std::vector<std::unique_ptr<TH1D>> universes;
    universes.reserve(nuniv);
    for (int k = 0; k < nuniv; ++k) {
        universes.emplace_back(
            make_total_mc_hist_weight_universe(spec, mc, weights_branch, k, "_univ"+std::to_string(k))
        );
    }
    return sample_covariance(*H0, universes);
}

TMatrixDSym cov_from_detvar_pm(const H1Spec& spec,
                               const std::vector<const Entry*>& mc,
                               const std::string& tag_up,
                               const std::string& tag_down) {
    auto H0   = make_total_mc_hist(spec, mc, "_nom");
    auto Hup  = make_total_mc_hist_detvar(spec, mc, tag_up,   "_up");
    auto Hdown= make_total_mc_hist_detvar(spec, mc, tag_down, "_down");
    return hessian_covariance(*H0, *Hup, *Hdown);
}

}
