#include "rarexsec/proc/Reweighter.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace rarexsec::proc {

Reweighter::Reweighter(std::vector<double> bin_edges,
                       std::vector<double> weights,
                       double default_weight)
    : bin_edges_(std::move(bin_edges)),
      weights_(std::move(weights)),
      default_weight_(default_weight) {
    if (bin_edges_.size() < 2) {
        throw std::invalid_argument("reweighter requires at least two bin edges");
    }
    if (bin_edges_.size() != weights_.size() + 1) {
        throw std::invalid_argument("reweighter bin edges and weights size mismatch");
    }
    if (!std::is_sorted(bin_edges_.begin(), bin_edges_.end())) {
        throw std::invalid_argument("reweighter bin edges must be sorted");
    }
}

double Reweighter::weight(double value) const {
    auto it = std::upper_bound(bin_edges_.begin(), bin_edges_.end(), value);
    if (it == bin_edges_.begin() || it == bin_edges_.end()) {
        return default_weight_;
    }
    std::size_t index = static_cast<std::size_t>(std::distance(bin_edges_.begin(), it) - 1);
    if (index >= weights_.size()) {
        return default_weight_;
    }
    double w = weights_[index];
    if (!std::isfinite(w) || w <= 0.0) {
        return default_weight_;
    }
    return w;
}

ROOT::RDF::RNode Reweighter::apply(ROOT::RDF::RNode node,
                                   const std::string& value_branch,
                                   const std::string& weight_branch) const {
    auto bin_edges = bin_edges_;
    auto weights = weights_;
    double default_weight = default_weight_;
    return node.Define(
        weight_branch,
        [bin_edges, weights, default_weight](double value) {
            auto it = std::upper_bound(bin_edges.begin(), bin_edges.end(), value);
            if (it == bin_edges.begin() || it == bin_edges.end()) {
                return default_weight;
            }
            std::size_t index = static_cast<std::size_t>(std::distance(bin_edges.begin(), it) - 1);
            if (index >= weights.size()) {
                return default_weight;
            }
            double w = weights[index];
            if (!std::isfinite(w) || w <= 0.0) {
                return default_weight;
            }
            return w;
        },
        {value_branch});
}

}
