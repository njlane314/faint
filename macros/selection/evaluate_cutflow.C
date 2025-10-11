#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>

#include <rarexsec/Hub.hh>
#include <rarexsec/proc/DataModel.hh>
#include <rarexsec/proc/Env.hh>
#include <rarexsec/proc/Selection.hh>

void evaluate_cutflow() {
    const auto env = rarexsec::Env::from_env();
    auto hub = env.make_hub();
    const auto mc = hub.simulation_entries(env.beamline, env.periods);

    auto is_signal = [](int ch_int) {
        const auto ch = static_cast<rarexsec::Channel>(ch_int);
        switch (ch) {
        case rarexsec::Channel::MuCC0pi_ge1p:
        case rarexsec::Channel::MuCC1pi:
        case rarexsec::Channel::MuCCPi0OrGamma:
        case rarexsec::Channel::MuCCNpi:
        case rarexsec::Channel::MuCCOther:
            return true;
        default:
            return false;
        }
    };

    auto sumw = [](ROOT::RDF::RNode n) {
        return static_cast<double>(n.Sum<float>("w_nominal").GetValue());
    };

    double denom = 0.0;
    for (const auto* rec : mc) {
        auto base = rec->nominal.rnode();
        denom += sumw(base.Filter([&](int ch){ return is_signal(ch); }, {"analysis_channels"}));
    }

    using Preset = rarexsec::selection::Preset;
    const std::vector<std::pair<std::string, Preset>> atoms = {
        {"Trigger",  Preset::Trigger},
        {"Slice",    Preset::Slice},
        {"Fiducial", Preset::Fiducial},
        {"Topology", Preset::Topology},
        {"Muon",     Preset::Muon}
    };

    std::string label;
    std::cout.setf(std::ios::fixed);
    std::cout.precision(6);
    std::cout << "Stage, Denom(signal), Selected(all), Selected(signal), Efficiency, Purity\n";

    for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (!label.empty()) label += "+";
        label += atoms[i].first;

        double sel_all = 0.0;
        double sel_sig = 0.0;

        for (const auto* rec : mc) {
            auto node = rec->nominal.rnode();
            for (std::size_t j = 0; j <= i; ++j)
                node = rarexsec::selection::apply(node, atoms[j].second, *rec);

            sel_all += sumw(node);
            sel_sig += sumw(node.Filter([&](int ch){ return is_signal(ch); }, {"analysis_channels"}));
        }

        const double eff = denom > 0.0 ? sel_sig / denom : 0.0;
        const double pur = sel_all > 0.0 ? sel_sig / sel_all : 0.0;

        std::cout << label << ", "
                  << denom << ", "
                  << sel_all << ", "
                  << sel_sig << ", "
                  << eff << ", "
                  << pur << "\n";
    }
}
