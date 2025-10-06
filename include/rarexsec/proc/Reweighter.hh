#pragma once
#include <ROOT/RDataFrame.hxx>
#include <string>
#include <vector>

namespace rarexsec::proc {

class Reweighter {
public:
    Reweighter(std::vector<double> bin_edges,
               std::vector<double> weights,
               double default_weight = 1.0);

    double weight(double value) const;

    ROOT::RDF::RNode apply(ROOT::RDF::RNode node,
                           const std::string& value_branch,
                           const std::string& weight_branch) const;

private:
    std::vector<double> bin_edges_;
    std::vector<double> weights_;
    double default_weight_;
};

}
