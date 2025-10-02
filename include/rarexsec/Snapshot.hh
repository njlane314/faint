#pragma once
#include "rarexsec/Hub.hh"

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace rarexsec {

namespace snapshot {

struct Options {
    std::string outdir = "snapshots";
    std::string tree = "analysis";
    std::vector<std::string> columns;
};

inline std::string origin_to_string(sample::origin k) {
    switch (k) {
        case sample::origin::data:
            return "data";
        case sample::origin::beam:
            return "beam";
        case sample::origin::strangeness:
            return "strangeness";
        case sample::origin::ext:
            return "ext";
        case sample::origin::dirt:
            return "dirt";
        default:
            return "unknown";
    }
}

inline std::string sanitise(std::string s) {
    for (char& c : s) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '_' || c == '.'))
            c = '_';
    }
    return s;
}

inline const std::vector<std::string>& default_columns() {
    static const std::vector<std::string> cols{
        "run", "subrun", "event",
        "w_nominal",
        "analysis_channels",
    };
    return cols;
}

inline std::vector<std::string> intersect_cols(ROOT::RDF::RNode node,
                                               const std::vector<std::string>& wanted) {
    auto have = node.GetColumnNames();
    std::unordered_set<std::string> avail(have.begin(), have.end());
    const auto& req = wanted.empty() ? default_columns() : wanted;

    std::vector<std::string> out;
    out.reserve(req.size());
    for (const auto& c : req) {
        if (avail.count(c))
            out.push_back(c);
    }
    return out;
}

inline std::string make_out_path(const Options& opt,
                                 const Entry& e,
                                 const std::string& detvar_tag) {
    const auto base = std::filesystem::path(e.file).filename().string();
    std::string name = sanitise(e.beamline) + "_" +
                       sanitise(e.period) + "_" +
                       sanitise(origin_to_string(e.kind));
    if (!detvar_tag.empty()) {
        name += "__" + sanitise(detvar_tag);
    }
    name += "__" + sanitise(base);
    name += ".root";

    std::filesystem::create_directories(opt.outdir);
    return (std::filesystem::path(opt.outdir) / name).string();
}

inline std::vector<std::string> write(const std::vector<const Entry*>& samples,
                                      const Options& opt = {}) {
    std::vector<std::string> outputs;
    outputs.reserve(samples.size());

    for (const Entry* e : samples) {
        if (!e)
            continue;

        const auto cols = intersect_cols(e->rnode(), opt.columns);
        const auto out = make_out_path(opt, *e, /*detvar*/ "");
        e->rnode().Snapshot(opt.tree, out, cols).GetValue();
        outputs.push_back(out);

        for (const auto& kv : e->detvars) {
            const auto& tag = kv.first;
            const auto& dv = kv.second;
            if (!dv.node)
                continue;

            const auto cols = intersect_cols(dv.rnode(), opt.columns);
            const auto out = make_out_path(opt, *e, tag);
            dv.rnode().Snapshot(opt.tree, out, cols).GetValue();
            outputs.push_back(out);
        }
    }

    return outputs;
}

std::vector<std::string> write(const Hub& hub,
                               std::string_view beamline,
                               const std::vector<std::string>& periods,
                               const Options& opt = {});

}
}
