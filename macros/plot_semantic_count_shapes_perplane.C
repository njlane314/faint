// plot_semantic_count_shapes_perplane.C
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TPad.h>
#include <TColor.h>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <iostream>
#include <cmath>

#include "rarexsec/Hub.hh"
#include "rarexsec/proc/Selection.hh"

// Load user libs (unchanged)
static void load_libs(const char* extra_libs) {
  gSystem->Load("librarexsec");
  if (!extra_libs) return;
  std::string s(extra_libs);
  for (char& c : s) if (c==';' || c==' ') c=',';
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) if (!tok.empty()) {
    if (gSystem->Load(tok.c_str()) < 0) std::cerr << "Warning: failed to load '" << tok << "'\n";
  }
}

static double sem_count(const std::vector<int>& counts, int lab) {
  if (counts.empty() || lab < 0 || lab >= (int)counts.size()) return 0.0;
  return static_cast<double>(counts[lab]);
}
static void normalize_pdf(TH1D& h) { double A = h.Integral("width"); if (A>0.0) h.Scale(1.0/A); }

static int discover_num_labels(const std::vector<const rarexsec::Entry*>& samples, const std::string& col) {
  using namespace rarexsec;
  int nlabels = 0;
  for (const auto* e : samples) if (e) {
    auto n0 = selection::apply(e->rnode(), selection::Preset::Empty, *e);
    auto n  = n0.Define("_rx_nlab_", [](const std::vector<int>& v){ return (int)v.size(); }, {col});
    auto v  = n.Take<int>("_rx_nlab_").GetValue();
    for (int c : v) nlabels = std::max(nlabels, c);
    if (nlabels > 0) break;
  }
  return nlabels;
}

// ─────────────────────────────────────────────────────────────────────────────
// Requested legend labels and palette (15 classes). Index = semantic label id.
constexpr std::size_t kPaletteSize = 15;
const std::array<const char*, kPaletteSize> kLegendLabels = {
    "#emptyset",
    "Cosmic",      "#mu",        "e^{-}",     "#gamma",
    "#pi^{#pm}",   "#pi^{0}",    "n",         "p",
    "K^{#pm}",     "K^{0}",      "#Lambda",   "#Sigma^{#pm}",
    "#Sigma^{0}",  "Other"
};
const std::array<int, kPaletteSize> kPalette = {
    TColor::GetColor("#333333"),
    TColor::GetColor("#666666"),
    TColor::GetColor("#e41a1c"),
    TColor::GetColor("#377eb8"),
    TColor::GetColor("#4daf4a"),
    TColor::GetColor("#ff7f00"),
    TColor::GetColor("#984ea3"),
    TColor::GetColor("#ffff33"),
    TColor::GetColor("#1b9e77"),
    TColor::GetColor("#f781bf"),
    TColor::GetColor("#a65628"),
    TColor::GetColor("#66a61e"),
    TColor::GetColor("#e6ab02"),
    TColor::GetColor("#a6cee3"),
    TColor::GetColor("#b15928")
};
inline std::size_t safe_idx(int lab) {
  if (lab < 0) return 0;
  return static_cast<std::size_t>(lab) < kPaletteSize ? static_cast<std::size_t>(lab) : kPaletteSize - 1;
}
inline std::string label_text(int lab) {
  return kLegendLabels[safe_idx(lab)];
}
inline int color_for_label(int lab) {
  return kPalette[safe_idx(lab)];
}

// Log-spaced bin edges helper (for x in [xmin, xmax])
inline std::vector<double> make_log_edges(double xmin, double xmax, int bins_per_decade = 40) {
  const double lx = std::log10(xmin), ux = std::log10(xmax);
  const int nbins = std::max(1, static_cast<int>(std::round((ux - lx) * bins_per_decade)));
  std::vector<double> edges; edges.reserve(nbins + 1);
  for (int i = 0; i <= nbins; ++i)
    edges.push_back(std::pow(10.0, lx + (ux - lx) * (double(i) / nbins)));
  return edges;
}

void plot_semantic_count_shapes_perplane(const char* extra_libs = "",
                                         bool use_event_counts   = true,
                                         bool normalize_to_pdf   = true,
                                         int nbins               = 60,     // (unused: kept for API compatibility)
                                         double xmin             = 0.0,     // (unused)
                                         double xmax             = 80000.0, // (unused)
                                         const char* beamline    = "numi-fhc",
                                         const char* period      = "run1",
                                         const char* config_path = "data/samples.json") {
  ROOT::EnableImplicitMT();
  load_libs(extra_libs);

  // Global style (no stats box, ticks on both sides, no exponent by default)
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TGaxis::SetMaxDigits(4);

  rarexsec::Hub hub(config_path);
  const std::vector<std::string> periods = {period};
  const auto mc = hub.simulation_entries(beamline, periods);
  if (mc.empty()) throw std::runtime_error("No MC samples found");

  const std::string out_dir = std::string("plots/") + beamline + "/" + period + "/semantic_count_shapes_perplane";
  gSystem->mkdir(out_dir.c_str(), true);

  const char* Ucnt = use_event_counts ? "event_semantic_counts_u" : "slice_semantic_counts_u";
  const char* Vcnt = use_event_counts ? "event_semantic_counts_v" : "slice_semantic_counts_v";
  const char* Ycnt = use_event_counts ? "event_semantic_counts_w" : "slice_semantic_counts_w"; // 'w' branch, but plane label is Y

  struct Plane { const char* tag; const char* col; const char* pretty; };
  const Plane planes[3] = { {"u", Ucnt, "U"}, {"v", Vcnt, "V"}, {"y", Ycnt, "Y"} };

  const int nlabels = discover_num_labels(mc, Ucnt);
  if (nlabels <= 0) { std::cerr << "No labels found.\n"; return; }

  for (const auto& p : planes) {
    std::vector<std::unique_ptr<TH1D>> H(nlabels);
    // Finer, log-spaced binning in [1, 1e4]
    constexpr int kBinsPerDecade = 40;
    const std::vector<double> log_edges = make_log_edges(1.0, 1e4, kBinsPerDecade);

    // pool all MC entries, no channel split
    for (size_t ie = 0; ie < mc.size(); ++ie) {
      const rarexsec::Entry* e = mc[ie];
      if (!e) continue;
      auto n0 = rarexsec::selection::apply(e->rnode(), rarexsec::selection::Preset::Empty, *e);

      for (int lab = 0; lab < nlabels; ++lab) {
        const std::string col = std::string("_rx_sc_") + p.tag + "_" + std::to_string(lab);
        auto n1 = n0.Define(col, [lab](const std::vector<int>& v){ return sem_count(v, lab); }, {p.col});

        ROOT::RDF::TH1DModel model(
            ("h_"+std::string(p.tag)+"_lab"+std::to_string(lab)+"_src"+std::to_string(ie)).c_str(),
            (";"+std::string(p.pretty)+"-plane semantic COUNT;Events").c_str(),
            log_edges);
        auto rr = n1.Histo1D(model, col, "w_nominal");

        const TH1D& part = rr.GetValue();
        if (!H[lab]) {
          H[lab].reset(static_cast<TH1D*>(part.Clone(("h_"+std::string(p.tag)+"_lab"+std::to_string(lab)).c_str())));
          H[lab]->SetDirectory(nullptr);
        } else {
          H[lab]->Add(&part);
        }
      }
    }

    // Style & (optional) PDF normalisation
    double ymax = 0.0;
    for (int lab = 0; lab < nlabels; ++lab) {
      if (!H[lab]) continue;
      if (normalize_to_pdf) normalize_pdf(*H[lab]);
      H[lab]->SetStats(false);
      H[lab]->SetFillStyle(0);
      H[lab]->SetLineColor(color_for_label(lab));
      H[lab]->SetLineWidth(3);
      ymax = std::max(ymax, H[lab]->GetMaximum());
    }

    const std::string cname = std::string("c_cnt_") + p.tag;
    TCanvas c(cname.c_str(), cname.c_str(), 900, 700);
    // Split canvas: legend on top, plot below.
    TPad* pad_main = new TPad(("pad_main_"+std::string(p.tag)).c_str(),
                              "pad_main", 0., 0.00, 1., 0.78);
    TPad* pad_leg  = new TPad(("pad_leg_"+std::string(p.tag)).c_str(),
                              "pad_leg",  0., 0.78, 1., 1.00);
    pad_main->SetLeftMargin(0.14);
    pad_main->SetRightMargin(0.05);
    pad_main->SetTopMargin(0.02);
    pad_main->SetBottomMargin(0.12);
    pad_main->SetLogx();               // requested log-x
    pad_leg->SetTopMargin(0.06);
    pad_leg->SetBottomMargin(0.10);
    pad_leg->SetLeftMargin(0.02);
    pad_leg->SetRightMargin(0.02);
    pad_main->Draw(); pad_leg->Draw();
    pad_main->cd();

    TH1D* frame = nullptr;
    for (int lab = 0; lab < nlabels; ++lab) if (H[lab]) { frame = H[lab].get(); break; }
    if (!frame) continue;

    // No plot title — only axis titles. Also: no exponents on axes.
    frame->SetTitle((";"+std::string(p.pretty)+"-plane semantic COUNT;"+std::string(normalize_to_pdf ? "Probability density" : "Events")).c_str());
    frame->GetXaxis()->SetNoExponent(true);
    frame->GetXaxis()->SetMoreLogLabels(true);
    frame->GetYaxis()->SetNoExponent(true);
    frame->GetXaxis()->SetMaxDigits(4);
    frame->GetYaxis()->SetMaxDigits(4);

    frame->SetMaximum( (ymax>0 ? 1.35*ymax : 1.) );
    frame->SetMinimum(0.0);
    // Show exactly [1, 1e4] on the x-axis
    frame->GetXaxis()->SetRangeUser(1., 1e4);
    frame->Draw("HIST");

    for (int lab = 0; lab < nlabels; ++lab)
      if (H[lab].get() != frame) H[lab]->Draw("HIST SAME");

    // Legend in its own pad (top), using requested labels & colors
    pad_leg->cd();
    TLegend leg(0.02, 0.05, 0.98, 0.95);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.035);
    leg.SetNColumns(nlabels > 10 ? 4 : (nlabels > 6 ? 3 : 2));
    for (int lab = 0; lab < nlabels; ++lab) if (H[lab]) {
      // Ensure legend line style/color match the plot
      H[lab]->SetLineColor(color_for_label(lab));
      leg.AddEntry(H[lab].get(), label_text(lab).c_str(), "l");
    }
    leg.Draw();
    pad_main->cd();

    // No TLatex banner/title here.

    c.SaveAs((out_dir + "/semantic_counts_" + p.tag + ".pdf").c_str());
  }

  std::cout << "Saved all-label semantic COUNT plots (no channel split) in: " << out_dir << std::endl;
}
