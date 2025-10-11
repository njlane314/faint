#pragma once
#include "rarexsec/Hub.hh"

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>

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
    std::string outfile = "all_samples.root";
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

inline std::string make_tree_name(const Options& opt, const Entry& e, const std::string& detvar) {
    std::string name = sanitise(opt.tree) + "__" +
                       sanitise(e.beamline) + "_" +
                       sanitise(e.period) + "_" +
                       sanitise(origin_to_string(e.kind));
    if (!detvar.empty())
        name += "__" + sanitise(detvar);
    return name;
}

inline std::string make_out_file(const Options& opt) {
    std::filesystem::create_directories(opt.outdir);
    return (std::filesystem::path(opt.outdir) / opt.outfile).string();
}

inline std::vector<std::string> write(const std::vector<const Entry*>& samples,
                                      const Options& opt = {}) {
    std::vector<std::string> outputs;
    outputs.reserve(1);

    const std::string outFile = make_out_file(opt);
    bool fileExists = std::filesystem::exists(outFile);

    auto snapshot_once = [&](ROOT::RDF::RNode node,
                             const std::string& treeName,
                             const std::vector<std::string>& cols) {
        ROOT::RDF::RSnapshotOptions sopt;
        sopt.fMode = fileExists ? "UPDATE" : "RECREATE";
        sopt.fOverwriteIfExists = true;
        node.Snapshot(treeName, outFile, cols, sopt).GetValue();
        fileExists = true;
    };

    for (const Entry* e : samples) {
        if (!e)
            continue;

        {
            const auto cols = intersect_cols(e->rnode(), opt.columns);
            const auto treeName = make_tree_name(opt, *e, "");
            snapshot_once(e->rnode(), treeName, cols);
        }

        for (const auto& kv : e->detvars) {
            const auto& tag = kv.first;
            const auto& dv = kv.second;
            if (!dv.node)
                continue;

            const auto cols = intersect_cols(dv.rnode(), opt.columns);
            const auto treeName = make_tree_name(opt, *e, tag);
            snapshot_once(dv.rnode(), treeName, cols);
        }
    }

    if (!samples.empty())
        outputs.push_back(outFile);
    return outputs;
}
}
}
