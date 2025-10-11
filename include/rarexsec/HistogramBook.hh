#pragma once

#include "rarexsec/DataModel.hh"
#include <TH1D.h>
#include <string>
#include <string_view>
#include <unordered_map>

namespace rarexsec::book {

inline ROOT::RDF::RResultPtr<TH1D>
H1(const Frame& f, const TH1D& model, std::string_view col,
   std::string_view wcol = "w_nominal") {
   return f.node.Histo1D(model, std::string(col), std::string(wcol));
}

inline std::unordered_map<std::string, ROOT::RDF::RResultPtr<TH1D>>
H1_variations(const Entry& e, const TH1D& model, std::string_view col,
              std::string_view wcol = "w_nominal") {
    std::unordered_map<std::string, ROOT::RDF::RResultPtr<TH1D>> out;
    out.reserve(e.detvars.size());
    for (auto const& kv : e.detvars) {
        out.emplace(kv.first, H1(kv.second, model, col, wcol));
    }
    return out;
}

}
