#pragma once
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>

namespace rarexsec {

namespace sample {
enum class origin { data, beam, strangeness, ext, dirt, unknown };

inline origin origin_from(const std::string& s) {
    if (s == "data")       return origin::data;
    if (s == "ext")          return origin::ext;
    if (s == "dirt")         return origin::dirt;
    if (s == "beam")         return origin::beam;
    if (s == "strangeness")  return origin::strangeness;
    return origin::unknown;
}
}

struct Data {
    std::shared_ptr<ROOT::RDataFrame> df;
    std::optional<ROOT::RDF::RNode> node;

    ROOT::RDF::RNode rnode() const { return ROOT::RDF::RNode{node.value()}; }
};

struct Entry {
    std::string beamline;
    std::string period;
    sample::origin kind = sample::origin::unknown;
    std::string file;
    double pot_nom  = 0.0;
    double pot_eqv  = 0.0;
    double trig_nom = 0.0;
    double trig_eqv = 0.0;
    Data nominal;
    std::unordered_map<std::string, Data> detvars;
    ROOT::RDF::RNode rnode() const { return nominal.rnode(); }
    const Data* detvar(const std::string& tag) const {
        auto it = detvars.find(tag);
        return it == detvars.end() ? nullptr : &it->second;
    }
};

class Hub {
public:
    explicit Hub(const std::string& path);

    std::vector<const rarexsec::Entry*> simulation(const std::string& beamline,
                                                   const std::vector<std::string>& periods) const;


private:
    using period_map   = std::unordered_map<std::string, std::vector<rarexsec::Entry>>;
    using beamline_map = std::unordered_map<std::string, period_map>;
    beamline_map db_;

    static Data sample(const Entry& rec);
};

}
