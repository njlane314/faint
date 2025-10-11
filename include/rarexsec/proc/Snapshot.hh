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

inline std::string source_to_string(Source s) {
    switch (s) {
        case Source::Data:
            return "data";
        case Source::Ext:
            return "ext";
        case Source::MC:
            return "mc";
    }
    return "unknown";
}

inline std::string slice_to_string(Slice s) {
    switch (s) {
        case Slice::None:
            return "none";
        case Slice::BeamInclusive:
            return "beam";
        case Slice::StrangenessInclusive:
            return "strangeness";
    }
    return "unknown";
}

inline std::string sample_label(const Entry& e) {
    if (e.kind == sample::origin::dirt)
        return "dirt";
    if (e.source == Source::MC) {
        const auto slice = slice_to_string(e.slice);
        if (slice == "none") return "mc";
        return slice;
    }
    return source_to_string(e.source);
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

inline std::vector<std::string> intersect_cols(ROOT::RDF::RNode node, const std::vector<std::string>& wanted) {
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

inline std::string make_out_path(const Options& opt, const Entry& e, const std::string& detvar) {
    const auto base = e.files.empty() ? std::string{}
                                      : std::filesystem::path(e.files.front()).filename().string();
    std::string name = sanitise(e.beamline) + "_" +
                       sanitise(e.period) + "_" +
                       sanitise(sample_label(e));
    if (!detvar.empty()) {
        name += "__" + sanitise(detvar);
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
            const auto cols = intersect_cols(dv.rnode(), opt.columns);
            const auto out = make_out_path(opt, *e, tag);
            dv.rnode().Snapshot(opt.tree, out, cols).GetValue();
            outputs.push_back(out);
        }
    }

    return outputs;
}

}
}