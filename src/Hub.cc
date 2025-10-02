#include "rarexsec/Hub.hh"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <fstream>
#include <stdexcept>

#include "rarexsec/Processor.hh"

rarexsec::Data rarexsec::Hub::sample(const Entry& rec) {
    constexpr const char* tree = "nuselection/EventSelectionFilter";

    auto df_ptr = std::make_shared<ROOT::RDataFrame>(tree, rec.file);
    ROOT::RDF::RNode node = *df_ptr;

    node = processor().run(node, rec);

    if (rec.kind == sample::origin::beam) node = node.Filter("!is_strange");
    else if (rec.kind == sample::origin::strangeness) node = node.Filter("is_strange");

    return Data{df_ptr, node};
}

rarexsec::Hub::Hub(const std::string& path) {
    std::ifstream cfg(path);
    if (!cfg) throw std::runtime_error("cannot open " + path);
    nlohmann::json j; cfg >> j;

    const auto& bl = j.at("beamlines");
    for (auto it_bl = bl.begin(); it_bl != bl.end(); ++it_bl) {
        const std::string beamline = it_bl.key();
        const auto& runs = it_bl.value();

        for (auto it_r = runs.begin(); it_r != runs.end(); ++it_r) {
            const std::string period = it_r.key();
            const auto& arr = it_r.value().at("samples");
            auto& bucket = db_[beamline][period];

            for (const auto& s : arr) {
                Entry rec;
                rec.beamline = beamline;
                rec.period   = period;
                rec.kind     = sample::origin_from(s.at("kind").get<std::string>());
                rec.file     = s.at("file").get<std::string>();
                
                rec.pot      = s.value("pot", 0.0);
                rec.pot_eff  = s.value("pot_eff", 0.0);
                rec.trig     = s.value("trig", 0.0);
                rec.trig_eff = s.value("trig_eff", 0.0);

                rec.nominal = sample(rec);

                if (s.contains("detvars")) {
                    const auto& dvs = s.at("detvars");
                    for (auto it_dv = dvs.begin(); it_dv != dvs.end(); ++it_dv) {
                        const std::string tag = it_dv.key();
                        const auto& desc = it_dv.value();
                        const std::string dv_file = desc.at("file").get<std::string>();
                        if (!dv_file.empty()) {
                            rec.detvars.emplace(tag, sample(dv_file, rec.kind, rec));
                        }
                    }
                }

                bucket.push_back(std::move(rec));
            }
        }
    }
}

rarexsec::Data rarexsec::Hub::sample(const std::string& file, sample::origin kind, const Entry& prototype) {
    Entry rec = prototype;
    rec.file = file;
    rec.kind = kind;
    rec.nominal = {};
    rec.detvars.clear();
    return sample(rec);
}

std::vector<const rarexsec::Entry*> rarexsec::Hub::simulation(const std::string& beamline,
                                                              const std::vector<std::string>& periods) const {
    std::vector<const Entry*> out;
    auto it_bl = db_.find(beamline);
    if (it_bl == db_.end()) return out;
    for (const auto& per : periods) {
        auto it_p = it_bl->second.find(per);
        if (it_p == it_bl->second.end()) continue;
        for (const auto& rec : it_p->second) {
            if (rec.kind != sample::origin::data) out.push_back(&rec);
        }
    }
    return out;
}
