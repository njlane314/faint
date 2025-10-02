#pragma once
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>

namespace rarexsec {

namespace sample {
enum class origin { data, beam, strangeness, ext, dirt, unknown };
}

struct Data {
    std::shared_ptr<ROOT::RDataFrame> df;
    ROOT::RDF::RNode node;
};

struct Entry {
    std::string beamline;
    std::string period;
    sample::origin kind = sample::origin::unknown;
    std::string file;
    double pot = 0.0;
    double trig = 0.0;
    Data nominal;
    std::unordered_map<std::string, Data> detvars;
    const ROOT::RDF::RNode& rnode() const { return nominal.node; }
    const Data* detvar(const std::string& tag) const {
        auto it = detvars.find(tag);
        return it == detvars.end() ? nullptr : &it->second;
    }
};

class Hub {
public:
    explicit Hub(const std::string& path);

    std::vector<const sample::Entry*> entries(const std::string& beamline,
                                                const std::string& period) const;

    std::vector<const sample::Entry*> simulation(const std::string& beamline,
                                                const std::vector<std::string>& periods) const;

    std::vector<std::string> beamlines() const;
    std::vector<std::string> periods(const std::string& beamline) const;

private:
    using period_map   = std::unordered_map<std::string, std::vector<sample::Entry>>;
    using beamline_map = std::unordered_map<std::string, period_map>;
    beamline_map db_;

    static sample::origin origin_from(const std::string& s);
    static bool is_simulation(sample::origin k);

    static Data make_frame(const std::vector<std::string>& files,
                            sample::origin kind);
};

} // namespace rarexsec
