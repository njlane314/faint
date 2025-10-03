#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TH1F.h>
#include <THStack.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLine.h>
#include <TColor.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <rarexsec/Hub.hh>
#include <rarexsec/proc/Selection.hh>
#include <rarexsec/Plotter.hh>
#include <rarexsec/plot/Channels.hh>

namespace rx = rarexsec;

struct StageDesc { rx::selection::Preset preset; const char* name; };

static inline void rx_set_style() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);
    gStyle->SetLabelSize(0.045, "XY");
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(1.10, "Y");
}

static inline void rx_ensure_dir(const std::string& d) {
    gSystem->mkdir(d.c_str(), kTRUE);
}

static inline std::string rx_fmt(double x, int prec = 2) {
    std::ostringstream s;
    s.setf(std::ios::fixed);
    s << std::setprecision(prec) << x;
    return s.str();
}

static inline const std::vector<StageDesc>& rx_stages() {
    static const std::vector<StageDesc> v = {
        {rx::selection::Preset::Empty, "Empty"},
        {rx::selection::Preset::Trigger, "Trigger"},
        {rx::selection::Preset::Slice, "Slice"},
        {rx::selection::Preset::Fiducial, "Fiducial (Reco)"},
        {rx::selection::Preset::Topology, "Topology"},
        {rx::selection::Preset::Muon, "Muon ID"},
        {rx::selection::Preset::InclusiveMuCC, "Inclusive #mu#nu CC (cum.)"}
    };
    return v;
}

static inline bool is_ext(rx::sample::origin k) {
    return k == rx::sample::origin::ext;
}

static inline const char* origin_to_str(rx::sample::origin k) {
    switch (k) {
        case rx::sample::origin::data: return "data";
        case rx::sample::origin::beam: return "beam";
        case rx::sample::origin::strangeness: return "strangeness";
        case rx::sample::origin::ext: return "ext";
        case rx::sample::origin::dirt: return "dirt";
        default: return "unknown";
    }
}

static inline int origin_color(rx::sample::origin k) {
    switch (k) {
        case rx::sample::origin::beam: return kAzure-4;
        case rx::sample::origin::strangeness: return kTeal+2;
        case rx::sample::origin::dirt: return kGray+2;
        case rx::sample::origin::ext: return kOrange+7;
        default: return kViolet-5;
    }
}

void selection_efficiency_purity_plots(const std::string& config_path = "data/samples.json",
                                       const std::string& beamline = "numi-fhc",
                                       const std::vector<std::string>& periods = {"run1"},
                                       const std::string& out_dir = "plots/selection",
                                       const std::string& img_fmt = "png") {
    try {
        ROOT::EnableImplicitMT();
        if (gSystem->Load("librarexsec") < 0) throw std::runtime_error("Failed to load librexsec");
        rx_set_style();
        rx_ensure_dir(out_dir);
        rx::Hub hub(config_path);
        const auto mc = hub.simulation_entries(beamline, periods);
        const auto data = hub.data_entries(beamline, periods);
        if (mc.empty() && data.empty()) throw std::runtime_error("No samples found");
        const auto& stages = rx_stages();
        const int NS = static_cast<int>(stages.size());
        auto truth_inclusive = [](int pdg, int ccnc, bool infv) { return (std::abs(pdg) == 14) && (ccnc == 0) && infv; };
        double denom_inclusive = 0.0, denom_actual = 0.0;
        for (const auto* e : mc) {
            if (!e) continue;
            auto n0 = e->rnode();
            if (!is_ext(e->kind)) {
                auto nt = n0.Filter(truth_inclusive, {"neutrino_pdg", "interaction_ccnc", "in_fiducial"});
                denom_inclusive += static_cast<double>(nt.Sum("w_nominal").GetValue());
            }
            auto na = n0.Filter("is_signal");
            denom_actual += static_cast<double>(na.Sum("w_nominal").GetValue());
        }
        std::vector<double> mc_total(NS, 0.0), S_incl(NS, 0.0), S_act(NS, 0.0);
        std::vector<ULong64_t> data_cnt(NS, 0ULL);
        for (int i = 0; i < NS; ++i) {
            double sum_mc = 0.0, s_inc = 0.0, s_act = 0.0;
            ULong64_t sum_d = 0ULL;
            for (const auto* e : mc) {
                if (!e) continue;
                auto n = e->rnode();
                for (int j = 0; j <= i; ++j) n = rx::selection::apply(n, stages[j].preset, *e);
                sum_mc += static_cast<double>(n.Sum("w_nominal").GetValue());
                if (!is_ext(e->kind)) {
                    auto ns = n.Filter(truth_inclusive, {"neutrino_pdg", "interaction_ccnc", "in_fiducial"});
                    s_inc += static_cast<double>(ns.Sum("w_nominal").GetValue());
                }
                auto na = n.Filter("is_signal");
                s_act += static_cast<double>(na.Sum("w_nominal").GetValue());
            }
            for (const auto* e : data) {
                if (!e) continue;
                auto n = e->rnode();
                for (int j = 0; j <= i; ++j) n = rx::selection::apply(n, stages[j].preset, *e);
                sum_d += static_cast<ULong64_t>(n.Count().GetValue());
            }
            mc_total[i] = sum_mc;
            data_cnt[i] = sum_d;
            S_incl[i] = s_inc;
            S_act[i] = s_act;
        }
        std::vector<double> purity_incl(NS, 0.0), purity_act(NS, 0.0), eff_incl(NS, 0.0), eff_act(NS, 0.0);
        for (int i = 0; i < NS; ++i) {
            purity_incl[i] = (mc_total[i] > 0.0) ? 100.0 * S_incl[i] / mc_total[i] : 0.0;
            purity_act[i] = (mc_total[i] > 0.0) ? 100.0 * S_act[i] / mc_total[i] : 0.0;
            eff_incl[i] = (denom_inclusive > 0.0) ? 100.0 * S_incl[i] / denom_inclusive : 0.0;
            eff_act[i] = (denom_actual > 0.0) ? 100.0 * S_act[i] / denom_actual : 0.0;
        }
        auto make_frame = [&](const char* name, const char* ytitle, double ymin, double ymax) {
            auto* h = new TH1F(name, "", NS, -0.5, NS - 0.5);
            for (int i = 0; i < NS; ++i) h->GetXaxis()->SetBinLabel(i + 1, stages[i].name);
            h->GetXaxis()->LabelsOption("v");
            h->SetMinimum(ymin);
            h->SetMaximum(ymax);
            h->GetYaxis()->SetTitle(ytitle);
            h->GetXaxis()->SetTitle("Selection stage (cumulative)");
            return h;
        };
        auto make_graph = [&](const std::vector<double>& y, int mstyle, int lstyle) {
            std::vector<double> x(NS);
            for (int i = 0; i < NS; ++i) x[i] = i;
            auto* g = new TGraph(NS, x.data(), y.data());
            g->SetMarkerStyle(mstyle);
            g->SetMarkerSize(1.0);
            g->SetLineStyle(lstyle);
            g->SetLineWidth(2);
            return g;
        };
        {
            TCanvas c("c_eff", "Signal efficiency vs stage", 900, 600);
            auto* frame = make_frame("frame_eff", "Efficiency / Acceptance [%]", 0.0, 100.0);
            frame->Draw("axis");
            auto* g_incl = make_graph(eff_incl, 20, 1);
            auto* g_act = make_graph(eff_act, 24, 2);
            g_incl->Draw("LP SAME");
            g_act->Draw("LP SAME");
            TLegend leg(0.60, 0.78, 0.93, 0.93);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(g_incl, "Inclusive #nu_{#mu} CC truth", "lp");
            leg.AddEntry(g_act, "Actual signal (strangeness)", "lp");
            leg.Draw();
            c.SaveAs((out_dir + "/efficiency_inclusive_vs_actual." + img_fmt).c_str());
            c.SaveAs((out_dir + "/efficiency_inclusive_vs_actual.pdf").c_str());
        }
        {
            TCanvas c("c_pur", "Purity vs stage", 900, 600);
            auto* frame = make_frame("frame_pur", "Purity [%]", 0.0, 100.0);
            frame->Draw("axis");
            auto* g_pincl = make_graph(purity_incl, 20, 1);
            auto* g_pact = make_graph(purity_act, 24, 2);
            g_pincl->Draw("LP SAME");
            g_pact->Draw("LP SAME");
            TLegend leg(0.60, 0.78, 0.93, 0.93);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(g_pincl, "Inclusive #nu_{#mu} CC purity", "lp");
            leg.AddEntry(g_pact, "Actual-signal purity", "lp");
            leg.Draw();
            c.SaveAs((out_dir + "/purity_inclusive_vs_actual." + img_fmt).c_str());
            c.SaveAs((out_dir + "/purity_inclusive_vs_actual.pdf").c_str());
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error in selection_efficiency_purity_plots: " << ex.what() << std::endl;
    }
}

void composition_by_stage(const std::string& config_path = "data/samples.json",
                          const std::string& beamline = "numi-fhc",
                          const std::vector<std::string>& periods = {"run1"},
                          const std::string& out_dir = "plots/composition",
                          const std::string& img_fmt = "png") {
    try {
        ROOT::EnableImplicitMT();
        if (gSystem->Load("librarexsec") < 0) throw std::runtime_error("Failed to load librexsec");
        rx_set_style();
        rx_ensure_dir(out_dir);
        rx::Hub hub(config_path);
        const auto mc = hub.simulation_entries(beamline, periods);
        const auto& stages = rx_stages();
        const int NS = static_cast<int>(stages.size());
        std::vector<rx::sample::origin> origins = {
            rx::sample::origin::beam, rx::sample::origin::strangeness,
            rx::sample::origin::dirt, rx::sample::origin::ext
        };
        std::map<rx::sample::origin, std::vector<double>> yield_by_origin;
        for (auto o : origins) yield_by_origin[o] = std::vector<double>(NS, 0.0);
        for (int i = 0; i < NS; ++i) {
            for (const auto* e : mc) {
                if (!e) continue;
                auto n = e->rnode();
                for (int j = 0; j <= i; ++j) n = rx::selection::apply(n, stages[j].preset, *e);
                yield_by_origin[e->kind][i] += static_cast<double>(n.Sum("w_nominal").GetValue());
            }
        }
        {
            TCanvas c("c_comp_origin", "Background composition by origin vs stage", 1000, 650);
            THStack hs("hs_origin", "Background composition by origin;Selection stage;Weighted yield");
            std::map<rx::sample::origin, TH1F*> hmap;
            for (auto o : origins) {
                auto* h = new TH1F((std::string("h_") + origin_to_str(o)).c_str(), "", NS, -0.5, NS - 0.5);
                for (int i = 0; i < NS; ++i) h->GetXaxis()->SetBinLabel(i + 1, stages[i].name);
                for (int i = 0; i < NS; ++i) h->SetBinContent(i + 1, yield_by_origin[o][i]);
                h->SetFillColor(origin_color(o));
                h->SetLineColor(kBlack);
                h->SetLineWidth(1);
                hs.Add(h, "hist");
                hmap[o] = h;
            }
            hs.Draw("nostack");
            hs.Draw("hist");
            hs.GetXaxis()->LabelsOption("v");
            TLegend leg(0.70, 0.70, 0.93, 0.93);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            for (auto o : origins) leg.AddEntry(hmap[o], origin_to_str(o), "f");
            leg.Draw();
            c.SaveAs((out_dir + "/composition_by_origin." + img_fmt).c_str());
            c.SaveAs((out_dir + "/composition_by_origin.pdf").c_str());
        }
        {
            const auto channels = rarexsec::plot::Channels::mc_keys();
            std::map<int, std::vector<double>> yield_by_channel;
            for (int ch : channels) yield_by_channel[ch] = std::vector<double>(NS, 0.0);
            for (int i = 0; i < NS; ++i) {
                for (const auto* e : mc) {
                    if (!e) continue;
                    auto n = e->rnode();
                    for (int j = 0; j <= i; ++j) n = rx::selection::apply(n, stages[j].preset, *e);
                    auto nc = n.Filter([](int c){ return c>=0; }, {"analysis_channels"});
                    for (int ch : channels) {
                        auto nf = nc.Filter([ch](int c){ return c==ch; }, {"analysis_channels"});
                        yield_by_channel[ch][i] += static_cast<double>(nf.Sum("w_nominal").GetValue());
                    }
                }
            }
            TCanvas c("c_comp_chan", "Background composition by analysis channel vs stage", 1200, 700);
            THStack hs("hs_chan", "Background composition by analysis channel;Selection stage;Weighted yield");
            std::vector<TH1F*> hlist;
            for (int ch : channels) {
                auto* h = new TH1F((std::string("h_ch_") + std::to_string(ch)).c_str(), "", NS, -0.5, NS - 0.5);
                for (int i = 0; i < NS; ++i) h->GetXaxis()->SetBinLabel(i + 1, stages[i].name);
                for (int i = 0; i < NS; ++i) h->SetBinContent(i + 1, yield_by_channel[ch][i]);
                h->SetFillColor(rarexsec::plot::Channels::color(ch));
                h->SetLineColor(kBlack);
                h->SetLineWidth(1);
                hs.Add(h, "hist");
                hlist.push_back(h);
            }
            hs.Draw("hist");
            hs.GetXaxis()->LabelsOption("v");
            TLegend leg(0.60, 0.10, 0.95, 0.90);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            int nadded = 0;
            for (int ch : channels) {
                if (nadded >= 18) break;
                leg.AddEntry(hlist[nadded], rarexsec::plot::Channels::label(ch).c_str(), "f");
                ++nadded;
            }
            leg.Draw();
            c.SaveAs((out_dir + "/composition_by_channel." + img_fmt).c_str());
            c.SaveAs((out_dir + "/composition_by_channel.pdf").c_str());
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error in composition_by_stage: " << ex.what() << std::endl;
    }
}

void scan_threshold_roc(const std::string& config_path = "data/samples.json",
                        const std::string& beamline = "numi-fhc",
                        const std::vector<std::string>& periods = {"run1"},
                        const std::string& param = "muon_min_track_score",
                        double start = 0.2, double stop = 0.9, double step = 0.05,
                        const std::string& out_dir = "plots/scans",
                        const std::string& img_fmt = "png") {
    try {
        ROOT::EnableImplicitMT();
        if (gSystem->Load("librarexsec") < 0) throw std::runtime_error("Failed to load librexsec");
        rx_set_style();
        rx_ensure_dir(out_dir);
        rx::Hub hub(config_path);
        const auto mc = hub.simulation_entries(beamline, periods);
        if (mc.empty()) throw std::runtime_error("No MC found");
        auto truth_inclusive = [](int pdg, int ccnc, bool infv) { return (std::abs(pdg) == 14) && (ccnc == 0) && infv; };
        double denom_inclusive = 0.0, denom_actual = 0.0;
        for (const auto* e : mc) {
            if (!e) continue;
            auto n0 = e->rnode();
            if (!is_ext(e->kind)) {
                auto nt = n0.Filter(truth_inclusive, {"neutrino_pdg", "interaction_ccnc", "in_fiducial"});
                denom_inclusive += static_cast<double>(nt.Sum("w_nominal").GetValue());
            }
            auto na0 = n0.Filter("is_signal");
            denom_actual += static_cast<double>(na0.Sum("w_nominal").GetValue());
        }
        std::vector<double> xs_eff_incl, ys_pur_incl, xs_eff_act, ys_pur_act, thr_vals;
        const bool is_muon_param = (param == "muon_min_track_score" || param == "muon_min_llr" || param == "muon_min_track_length" || param == "muon_max_track_distance");
        const bool is_topo_param = (param == "topology_min_contained_fraction" || param == "topology_min_cluster_fraction");
        if (!is_muon_param && !is_topo_param) throw std::runtime_error("Unsupported parameter");
        for (double thr = start; (step > 0 ? thr <= stop : thr >= stop); thr += step) {
            double tot_mc = 0.0, s_inc = 0.0, s_act = 0.0;
            for (const auto* e : mc) {
                if (!e) continue;
                auto n = e->rnode();
                n = rx::selection::apply(n, rx::selection::Preset::Trigger, *e);
                n = rx::selection::apply(n, rx::selection::Preset::Slice, *e);
                n = rx::selection::apply(n, rx::selection::Preset::Fiducial, *e);
                if (is_topo_param) {
                    const double cmin = (param == "topology_min_contained_fraction") ? thr : rx::selection::topology_min_contained_fraction;
                    const double lmin = (param == "topology_min_cluster_fraction") ? thr : rx::selection::topology_min_cluster_fraction;
                    n = n.Filter([=](float cf, float cl){ return cf >= cmin && cl >= lmin; }, {"contained_fraction", "slice_cluster_fraction"});
                    n = rx::selection::apply(n, rx::selection::Preset::Muon, *e);
                } else {
                    n = rx::selection::apply(n, rx::selection::Preset::Topology, *e);
                    const double smin = (param == "muon_min_track_score") ? thr : rx::selection::muon_min_track_score;
                    const double lmin = (param == "muon_min_llr") ? thr : rx::selection::muon_min_llr;
                    const double lenmin = (param == "muon_min_track_length") ? thr : rx::selection::muon_min_track_length;
                    const double dmax = (param == "muon_max_track_distance") ? thr : rx::selection::muon_max_track_distance;
                    const unsigned genreq = rx::selection::muon_required_generation;
                    n = n.Filter([=](const ROOT::RVec<float>& scores,
                                     const ROOT::RVec<float>& llrs,
                                     const ROOT::RVec<float>& lengths,
                                     const ROOT::RVec<float>& distances,
                                     const ROOT::RVec<unsigned>& generations) {
                        const auto n = scores.size();
                        for (std::size_t i = 0; i < n; ++i) {
                            const bool pass = scores[i] > smin && llrs[i] > lmin && lengths[i] > lenmin && distances[i] < dmax && generations[i] == genreq;
                            if (pass) return true;
                        }
                        return false;
                    }, {"track_shower_scores", "trk_llr_pid_v", "track_length", "track_distance_to_vertex", "pfp_generations"});
                }
                tot_mc += static_cast<double>(n.Sum("w_nominal").GetValue());
                if (!is_ext(e->kind)) {
                    auto ns = n.Filter(truth_inclusive, {"neutrino_pdg", "interaction_ccnc", "in_fiducial"});
                    s_inc += static_cast<double>(ns.Sum("w_nominal").GetValue());
                }
                auto na = n.Filter("is_signal");
                s_act += static_cast<double>(na.Sum("w_nominal").GetValue());
            }
            const double pur_incl = (tot_mc > 0.0) ? 100.0 * s_inc / tot_mc : 0.0;
            const double pur_act = (tot_mc > 0.0) ? 100.0 * s_act / tot_mc : 0.0;
            const double eff_incl = (denom_inclusive > 0.0) ? 100.0 * s_inc / denom_inclusive : 0.0;
            const double eff_act = (denom_actual > 0.0) ? 100.0 * s_act / denom_actual : 0.0;
            xs_eff_incl.push_back(eff_incl);
            ys_pur_incl.push_back(pur_incl);
            xs_eff_act.push_back(eff_act);
            ys_pur_act.push_back(pur_act);
            thr_vals.push_back(thr);
        }
        TCanvas c("c_roc", "ROC scan", 900, 650);
        auto* g_incl = new TGraph(xs_eff_incl.size(), xs_eff_incl.data(), ys_pur_incl.data());
        auto* g_act = new TGraph(xs_eff_act.size(), xs_eff_act.data(), ys_pur_act.data());
        auto* frame = new TH1F("frame", "", 10, 0, 100);
        frame->SetMinimum(0);
        frame->SetMaximum(100);
        frame->GetXaxis()->SetTitle("Efficiency [%]");
        frame->GetYaxis()->SetTitle("Purity [%]");
        frame->Draw("axis");
        g_incl->SetMarkerStyle(20);
        g_incl->SetLineStyle(1);
        g_incl->SetLineWidth(2);
        g_act->SetMarkerStyle(24);
        g_act->SetLineStyle(2);
        g_act->SetLineWidth(2);
        g_incl->Draw("LP SAME");
        g_act->Draw("LP SAME");
        TLegend leg(0.58, 0.78, 0.93, 0.93);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(g_incl, "Inclusive #nu_{#mu} CC truth", "lp");
        leg.AddEntry(g_act, "Actual signal (strangeness)", "lp");
        leg.Draw();
        TLatex ltx;
        ltx.SetTextSize(0.03);
        for (size_t i = 0; i < thr_vals.size(); ++i) {
            ltx.DrawLatex(xs_eff_incl[i] + 0.5, ys_pur_incl[i] + 0.5, (param + "=" + rx_fmt(thr_vals[i], 2)).c_str());
        }
        c.SaveAs((out_dir + "/roc_" + param + "." + img_fmt).c_str());
        c.SaveAs((out_dir + "/roc_" + param + ".pdf").c_str());
    } catch (const std::exception& ex) {
        std::cerr << "Error in scan_threshold_roc: " << ex.what() << std::endl;
    }
}

void evaluate_recognised_signal_vs_stage(const std::string& config_path = "data/samples.json",
                                         const std::string& beamline = "numi-fhc",
                                         const std::vector<std::string>& periods = {"run1"},
                                         const std::string& out_dir = "plots/recognised",
                                         const std::string& img_fmt = "png") {
    try {
        ROOT::EnableImplicitMT();
        if (gSystem->Load("librarexsec") < 0) throw std::runtime_error("Failed to load librexsec");
        rx_set_style();
        rx_ensure_dir(out_dir);
        rx::Hub hub(config_path);
        const auto mc = hub.simulation_entries(beamline, periods);
        const auto data = hub.data_entries(beamline, periods);
        if (mc.empty() && data.empty()) throw std::runtime_error("No samples found");
        const auto& stages = rx_stages();
        const int NS = static_cast<int>(stages.size());
        auto truth_inclusive = [](int pdg, int ccnc, bool infv) { return (std::abs(pdg) == 14) && (ccnc == 0) && infv; };
        double denom_inclusive = 0.0;
        for (const auto* e : mc) {
            if (!e) continue;
            if (!is_ext(e->kind)) {
                auto nt = e->rnode().Filter(truth_inclusive, {"neutrino_pdg", "interaction_ccnc", "in_fiducial"});
                denom_inclusive += static_cast<double>(nt.Sum("w_nominal").GetValue());
            }
        }
        std::vector<double> mc_total(NS, 0.0), S_rec(NS, 0.0);
        for (int i = 0; i < NS; ++i) {
            for (const auto* e : mc) {
                if (!e) continue;
                auto n = e->rnode();
                for (int j = 0; j <= i; ++j) n = rx::selection::apply(n, stages[j].preset, *e);
                mc_total[i] += static_cast<double>(n.Sum("w_nominal").GetValue());
                auto nr = n.Filter("recognised_signal");
                S_rec[i] += static_cast<double>(nr.Sum("w_nominal").GetValue());
            }
        }
        std::vector<double> purity_rec(NS, 0.0), eff_rec(NS, 0.0);
        for (int i = 0; i < NS; ++i) {
            purity_rec[i] = (mc_total[i] > 0.0) ? 100.0 * S_rec[i] / mc_total[i] : 0.0;
            eff_rec[i] = (denom_inclusive > 0.0) ? 100.0 * S_rec[i] / denom_inclusive : 0.0;
        }
        auto make_frame = [&](const char* name, const char* ytitle, double ymin, double ymax) {
            auto* h = new TH1F(name, "", NS, -0.5, NS - 0.5);
            for (int i = 0; i < NS; ++i) h->GetXaxis()->SetBinLabel(i + 1, stages[i].name);
            h->GetXaxis()->LabelsOption("v");
            h->SetMinimum(ymin);
            h->SetMaximum(ymax);
            h->GetYaxis()->SetTitle(ytitle);
            h->GetXaxis()->SetTitle("Selection stage (cumulative)");
            return h;
        };
        auto make_graph = [&](const std::vector<double>& y, int mstyle, int lstyle) {
            std::vector<double> x(NS);
            for (int i = 0; i < NS; ++i) x[i] = i;
            auto* g = new TGraph(NS, x.data(), y.data());
            g->SetMarkerStyle(mstyle);
            g->SetMarkerSize(1.0);
            g->SetLineStyle(lstyle);
            g->SetLineWidth(2);
            return g;
        };
        {
            TCanvas c("c_eff_rec", "Recognised signal efficiency", 900, 600);
            auto* frame = make_frame("frame_eff_rec", "Efficiency [%]", 0.0, 100.0);
            frame->Draw("axis");
            auto* g = make_graph(eff_rec, 20, 1);
            g->Draw("LP SAME");
            c.SaveAs((out_dir + "/recognised_efficiency." + img_fmt).c_str());
            c.SaveAs((out_dir + "/recognised_efficiency.pdf").c_str());
        }
        {
            TCanvas c("c_pur_rec", "Recognised signal purity", 900, 600);
            auto* frame = make_frame("frame_pur_rec", "Purity [%]", 0.0, 100.0);
            frame->Draw("axis");
            auto* g = make_graph(purity_rec, 24, 2);
            g->Draw("LP SAME");
            c.SaveAs((out_dir + "/recognised_purity." + img_fmt).c_str());
            c.SaveAs((out_dir + "/recognised_purity.pdf").c_str());
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error in evaluate_recognised_signal_vs_stage: " << ex.what() << std::endl;
    }
}

void plot_final_stage_kinematics(const std::string& config_path = "data/samples.json",
                                 const std::string& beamline = "numi-fhc",
                                 const std::vector<std::string>& periods = {"run1"},
                                 const std::string& out_dir = "plots/kinematics",
                                 const std::string& img_fmt = "png") {
    try {
        ROOT::EnableImplicitMT();
        if (gSystem->Load("librarexsec") < 0) throw std::runtime_error("Failed to load librexsec");
        rx_set_style();
        rx_ensure_dir(out_dir);
        rx::Hub hub(config_path);
        const auto mc = hub.simulation_entries(beamline, periods);
        const auto data = hub.data_entries(beamline, periods);
        rarexsec::plot::Options opt;
        opt.out_dir = out_dir;
        opt.image_format = img_fmt;
        opt.show_ratio = true;
        opt.legend_on_top = true;
        opt.legend_split = 0.82;
        opt.x_title = "";
        opt.y_title = "Events";
        opt.beamline = beamline;
        opt.periods = periods;
        rarexsec::plot::Plotter plotter(opt);
        std::vector<rarexsec::plot::H1Spec> hspecs;
        {
            rarexsec::plot::H1Spec s;
            s.id = "topological_score";
            s.title = ";Topological score;Events";
            s.expr = "topological_score";
            s.nbins = 20; s.xmin = 0; s.xmax = 1.0;
            s.sel = rx::selection::Preset::InclusiveMuCC;
            hspecs.push_back(s);
        }
        {
            rarexsec::plot::H1Spec s;
            s.id = "contained_fraction";
            s.title = ";Contained fraction;Events";
            s.expr = "contained_fraction";
            s.nbins = 20; s.xmin = 0; s.xmax = 1.0;
            s.sel = rx::selection::Preset::InclusiveMuCC;
            hspecs.push_back(s);
        }
        {
            rarexsec::plot::H1Spec s;
            s.id = "slice_cluster_fraction";
            s.title = ";Slice cluster fraction;Events";
            s.expr = "slice_cluster_fraction";
            s.nbins = 20; s.xmin = 0; s.xmax = 1.0;
            s.sel = rx::selection::Preset::InclusiveMuCC;
            hspecs.push_back(s);
        }
        {
            rarexsec::plot::H1Spec s;
            s.id = "reco_vtx_x";
            s.title = ";Reco vertex x [cm];Events";
            s.expr = "reco_neutrino_vertex_sce_x";
            s.nbins = 30; s.xmin = -200; s.xmax = 200;
            s.sel = rx::selection::Preset::InclusiveMuCC;
            hspecs.push_back(s);
        }
        {
            rarexsec::plot::H1Spec s;
            s.id = "reco_vtx_y";
            s.title = ";Reco vertex y [cm];Events";
            s.expr = "reco_neutrino_vertex_sce_y";
            s.nbins = 30; s.xmin = -200; s.xmax = 200;
            s.sel = rx::selection::Preset::InclusiveMuCC;
            hspecs.push_back(s);
        }
        {
            rarexsec::plot::H1Spec s;
            s.id = "reco_vtx_z";
            s.title = ";Reco vertex z [cm];Events";
            s.expr = "reco_neutrino_vertex_sce_z";
            s.nbins = 40; s.xmin = 0; s.xmax = 1100;
            s.sel = rx::selection::Preset::InclusiveMuCC;
            hspecs.push_back(s);
        }
        for (const auto& s : hspecs) plotter.draw_stack_by_channel(s, mc, data);
    } catch (const std::exception& ex) {
        std::cerr << "Error in plot_final_stage_kinematics: " << ex.what() << std::endl;
    }
}

void period_stability_summary(const std::string& config_path = "data/samples.json",
                              const std::string& beamline = "numi-fhc",
                              const std::vector<std::string>& periods = {"run1","run2","run3"},
                              const std::string& out_dir = "plots/stability",
                              const std::string& img_fmt = "png") {
    try {
        ROOT::EnableImplicitMT();
        if (gSystem->Load("librarexsec") < 0) throw std::runtime_error("Failed to load librexsec");
        rx_set_style();
        rx_ensure_dir(out_dir);
        rx::Hub hub(config_path);
        auto truth_inclusive = [](int pdg, int ccnc, bool infv) { return (std::abs(pdg) == 14) && (ccnc == 0) && infv; };
        std::vector<double> eff_incl, pur_incl, eff_act, pur_act;
        std::vector<std::string> labs;
        for (const auto& per : periods) {
            const auto mc = hub.simulation_entries(beamline, {per});
            const auto data = hub.data_entries(beamline, {per});
            double denom_inclusive = 0.0, denom_actual = 0.0;
            for (const auto* e : mc) {
                if (!e) continue;
                if (!is_ext(e->kind)) {
                    auto nt = e->rnode().Filter(truth_inclusive, {"neutrino_pdg", "interaction_ccnc", "in_fiducial"});
                    denom_inclusive += static_cast<double>(nt.Sum("w_nominal").GetValue());
                }
                auto na0 = e->rnode().Filter("is_signal");
                denom_actual += static_cast<double>(na0.Sum("w_nominal").GetValue());
            }
            double tot_mc = 0.0, s_inc = 0.0, s_act = 0.0;
            for (const auto* e : mc) {
                if (!e) continue;
                auto n = e->rnode();
                n = rx::selection::apply(n, rx::selection::Preset::InclusiveMuCC, *e);
                tot_mc += static_cast<double>(n.Sum("w_nominal").GetValue());
                if (!is_ext(e->kind)) {
                    auto ns = n.Filter(truth_inclusive, {"neutrino_pdg", "interaction_ccnc", "in_fiducial"});
                    s_inc += static_cast<double>(ns.Sum("w_nominal").GetValue());
                }
                auto na = n.Filter("is_signal");
                s_act += static_cast<double>(na.Sum("w_nominal").GetValue());
            }
            double p_incl = (tot_mc > 0.0) ? 100.0 * s_inc / tot_mc : 0.0;
            double p_act = (tot_mc > 0.0) ? 100.0 * s_act / tot_mc : 0.0;
            double e_incl = (denom_inclusive > 0.0) ? 100.0 * s_inc / denom_inclusive : 0.0;
            double e_act = (denom_actual > 0.0) ? 100.0 * s_act / denom_actual : 0.0;
            labs.push_back(per);
            pur_incl.push_back(p_incl);
            pur_act.push_back(p_act);
            eff_incl.push_back(e_incl);
            eff_act.push_back(e_act);
        }
        auto make_frame = [&](const char* name, const char* ytitle, int N, double ymin, double ymax) {
            auto* h = new TH1F(name, "", N, -0.5, N - 0.5);
            for (int i = 0; i < N; ++i) h->GetXaxis()->SetBinLabel(i + 1, labs[i].c_str());
            h->GetXaxis()->LabelsOption("v");
            h->SetMinimum(ymin);
            h->SetMaximum(ymax);
            h->GetYaxis()->SetTitle(ytitle);
            h->GetXaxis()->SetTitle("Period");
            return h;
        };
        auto make_graph = [&](const std::vector<double>& y, int mstyle, int lstyle) {
            const int N = static_cast<int>(y.size());
            std::vector<double> x(N);
            for (int i = 0; i < N; ++i) x[i] = i;
            auto* g = new TGraph(N, x.data(), y.data());
            g->SetMarkerStyle(mstyle);
            g->SetMarkerSize(1.0);
            g->SetLineStyle(lstyle);
            g->SetLineWidth(2);
            return g;
        };
        {
            const int N = static_cast<int>(labs.size());
            TCanvas c("c_stab_eff", "Period stability: efficiency", 950, 650);
            auto* frame = make_frame("frame_stab_eff", "Efficiency [%]", N, 0.0, 100.0);
            frame->Draw("axis");
            auto* g_incl = make_graph(eff_incl, 20, 1);
            auto* g_act = make_graph(eff_act, 24, 2);
            g_incl->Draw("LP SAME");
            g_act->Draw("LP SAME");
            TLegend leg(0.60, 0.78, 0.93, 0.93);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(g_incl, "Inclusive #nu_{#mu} CC", "lp");
            leg.AddEntry(g_act, "Actual signal", "lp");
            leg.Draw();
            c.SaveAs((out_dir + "/stability_efficiency." + img_fmt).c_str());
            c.SaveAs((out_dir + "/stability_efficiency.pdf").c_str());
        }
        {
            const int N = static_cast<int>(labs.size());
            TCanvas c("c_stab_pur", "Period stability: purity", 950, 650);
            auto* frame = make_frame("frame_stab_pur", "Purity [%]", N, 0.0, 100.0);
            frame->Draw("axis");
            auto* g_incl = make_graph(pur_incl, 20, 1);
            auto* g_act = make_graph(pur_act, 24, 2);
            g_incl->Draw("LP SAME");
            g_act->Draw("LP SAME");
            TLegend leg(0.60, 0.78, 0.93, 0.93);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(g_incl, "Inclusive #nu_{#mu} CC", "lp");
            leg.AddEntry(g_act, "Actual signal", "lp");
            leg.Draw();
            c.SaveAs((out_dir + "/stability_purity." + img_fmt).c_str());
            c.SaveAs((out_dir + "/stability_purity.pdf").c_str());
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error in period_stability_summary: " << ex.what() << std::endl;
    }
}
