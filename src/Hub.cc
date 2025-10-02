#include "Hub.hh"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <fstream>
#include <stdexcept>

using json = nlohmann::json;

static std::vector<std::string> vec_str(const json& v) {
    if (v.is_string()) return { v.get<std::string>() };
    if (v.is_array())  return v.get<std::vector<std::string>>();
    return {};
}

static std::vector<std::string> fetch_list(const json& obj, const char* key) {
    if (!obj.contains(key)) return {};
    return vec_str(obj.at(key));
}

rarexsec::sample::origin rarexsec::Hub::origin_from(const std::string& s) {
    if (s == "data")         return sample::origin::data;
    if (s == "ext")          return sample::origin::ext;
    if (s == "dirt")         return sample::origin::dirt;
    if (s == "beam")         return sample::origin::beam;
    if (s == "strangeness")  return sample::origin::strangeness;
    return sample::origin::unknown;
}

bool rarexsec::Hub::is_simulation(sample::origin k) {
    return k != sample::origin::data;
}

rarexsec::Frame rarexsec::Hub::make_frame(const std::vector<std::string>& files,
                                          sample::origin kind) {
    constexpr const char* kTree = "nuselection/EventSelectionFilter";
    constexpr const char* kTruthStrange = "is_strange";

    auto chain = std::make_shared<TChain>(kTree);
    for (const auto& f : files) chain->Add(f.c_str());

    auto df_ptr = std::make_shared<ROOT::RDataFrame>(*chain);
    ROOT::RDF::RNode node = *df_ptr;

    if (kind == sample::origin::beam) node = node.Filter(std::string(kTruthStrange) + " == 0", "truth_beam");
    else if (kind == sample::origin::strangeness) node = node.Filter(std::string(kTruthStrange) + " != 0", "truth_strangeness");

    return Frame{chain, df_ptr, node};
}

rarexsec::Hub::Hub(const std::string& path) {
    std::ifstream cfg(path);
    if (!cfg) throw std::runtime_error("cannot open " + path);
    json j; cfg >> j;

    const auto& bl = j.at("beamlines");
    for (auto it_bl = bl.begin(); it_bl != bl.end(); ++it_bl) {
        const std::string beamline = it_bl.key();
        const auto& runs = it_bl.value();

        for (auto it_r = runs.begin(); it_r != runs.end(); ++it_r) {
            const std::string period = it_r.key();
            const auto& arr = it_r.value().at("samples");
            auto& bucket = db_[beamline][period];

            for (const auto& s : arr) {
                sample::Entry rec;
                rec.beamline = beamline;
                rec.period   = period;
                rec.kind     = origin_from(s.at("kind").get<std::string>());
                rec.files    = s.at("files").get<std::vector<std::string>>();
                rec.pot      = s.value("pot", 0.0);

                rec.nominal = make_frame(rec.files, rec.kind);

                if (s.contains("detvars")) {
                    const auto& dvs = s.at("detvars");
                    for (auto it_dv = dvs.begin(); it_dv != dvs.end(); ++it_dv) {
                        const std::string tag = it_dv.key();
                        const auto& desc = it_dv.value();
                        std::vector<std::string> dv_files = fetch_list(desc, "files");
                        if (!dv_files.empty()) {
                        rec.detvars.emplace(tag, make_frame(dv_files, rec.kind));
                        }
                    }
                }

                bucket.push_back(std::move(rec));
            }
        }
    }
}

std::vector<const rarexsec::sample::Entry*> rarexsec::Hub::entries(const std::string& beamline,
                                                                   const std::string& period) const {
    std::vector<const sample::Entry*> out;
    auto it_bl = db_.find(beamline);
    if (it_bl == db_.end()) return out;
    auto it_p = it_bl->second.find(period);
    if (it_p == it_bl->second.end()) return out;
    for (const auto& rec : it_p->second) out.push_back(&rec);
    return out;
}

std::vector<const rarexsec::sample::Entry*> rarexsec::Hub::simulation(const std::string& beamline,
                                                                      const std::vector<std::string>& periods) const {
    std::vector<const sample::Entry*> out;
    auto it_bl = db_.find(beamline);
    if (it_bl == db_.end()) return out;
    for (const auto& per : periods) {
        auto it_p = it_bl->second.find(per);
        if (it_p == it_bl->second.end()) continue;
        for (const auto& rec : it_p->second) if (is_simulation(rec.kind)) out.push_back(&rec);
    }
    return out;
}

std::vector<std::string> rarexsec::Hub::beamlines() const {
    std::vector<std::string> v; v.reserve(db_.size());
    for (const auto& kv : db_) v.push_back(kv.first);
    std::sort(v.begin(), v.end());
    return v;
}

std::vector<std::string> rarexsec::Hub::periods(const std::string& beamline) const {
    std::vector<std::string> v;
    auto it = db_.find(beamline);
    if (it == db_.end()) return v;
    for (const auto& kv : it->second) v.push_back(kv.first);
    std::sort(v.begin(), v.end());
    return v;
}