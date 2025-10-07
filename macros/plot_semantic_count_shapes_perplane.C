// plot_semantic_count_shapes_perplane.C
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TLatex.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <iostream>

#include "rarexsec/Hub.hh"
#include "rarexsec/proc/Selection.hh"

// Load user libs
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
  for (const auto* e : samples) if (e) {
    auto n0 = selection::apply(e->rnode(), selection::Preset::Empty, *e);
    auto n  = n0.Define("_rx_nlab_", [](const std::vector<int>& v){ return (int)v.size(); }, {col});
    auto v  = n.Take<int>("_rx_nlab_").GetValue();
    if (!v.empty()) return v.front();
  }
  return 0;
}
// simple palette
static int color_for_label(int i){ static const int P[]={kBlack,kRed+1,kAzure+2,kGreen+2,kOrange+7,kMagenta+1,kCyan+2,kViolet+1,kTeal+3,kPink+7,kGray+2}; return P[i% (sizeof(P)/sizeof(int))]; }
static std::string label_name(int i){ std::ostringstream ss; ss<<"label "<<i; return ss.str(); }

void plot_semantic_count_shapes_perplane(const char* extra_libs = "",
                                         bool use_event_counts   = true,
                                         bool normalize_to_pdf   = true,
                                         int nbins               = 60,
                                         double xmin             = 0.0,
                                         double xmax             = 80000.0,
                                         const char* beamline    = "numi-fhc",
                                         const char* period      = "run1",
                                         const char* config_path = "data/samples.json") {
  ROOT::EnableImplicitMT();
  load_libs(extra_libs);

  rarexsec::Hub hub(config_path);
  const std::vector<std::string> periods = {period};
  const auto mc = hub.simulation_entries(beamline, periods);
  if (mc.empty()) throw std::runtime_error("No MC samples found");

  const std::string out_dir = std::string("plots/") + beamline + "/" + period + "/semantic_count_shapes_perplane";
  gSystem->mkdir(out_dir.c_str(), true);

  const char* Ucnt = use_event_counts ? "event_semantic_counts_u" : "slice_semantic_counts_u";
  const char* Vcnt = use_event_counts ? "event_semantic_counts_v" : "slice_semantic_counts_v";
  const char* Wcnt = use_event_counts ? "event_semantic_counts_w" : "slice_semantic_counts_w";

  struct Plane { const char* tag; const char* col; };
  const Plane planes[3] = { {"u", Ucnt}, {"v", Vcnt}, {"w", Wcnt} };

  const int nlabels = discover_num_labels(mc, Ucnt);
  if (nlabels <= 0) { std::cerr << "No labels found.\n"; return; }

  for (const auto& p : planes) {
    std::vector<std::unique_ptr<TH1D>> H(nlabels);
    // pool all MC entries, no channel split
    for (size_t ie = 0; ie < mc.size(); ++ie) {
      const rarexsec::Entry* e = mc[ie];
      if (!e) continue;
      auto n0 = rarexsec::selection::apply(e->rnode(), rarexsec::selection::Preset::Empty, *e);
      for (int lab = 0; lab < nlabels; ++lab) {
        const std::string col = std::string("_rx_sc_") + p.tag + "_" + std::to_string(lab);
    auto n1 = n0.Define(col, [lab](const std::vector<int>& v){ return sem_count(v, lab); }, {p.col});
    auto rr = n1.Histo1D(ROOT::RDF::TH1DModel(("h_"+std::string(p.tag)+"_lab"+std::to_string(lab)+"_src"+std::to_string(ie)).c_str(),
                                             (";"+std::string(1,(char)toupper(*p.tag))+"-plane semantic COUNT;Events").c_str(),
                                             nbins, xmin, xmax),
                          col, "w_nominal");
        const TH1D& part = rr.GetValue();
        if (!H[lab]) {
          H[lab].reset(static_cast<TH1D*>(part.Clone(("h_"+std::string(p.tag)+"_lab"+std::to_string(lab)).c_str())));
          H[lab]->SetDirectory(nullptr);
        } else H[lab]->Add(&part);
      }
    }

    double ymax = 0.0;
    for (int lab = 0; lab < nlabels; ++lab) {
      if (!H[lab]) continue;
      if (normalize_to_pdf) normalize_pdf(*H[lab]);
      H[lab]->SetFillStyle(0);
      H[lab]->SetLineColor(color_for_label(lab));
      H[lab]->SetLineWidth(3);
      ymax = std::max(ymax, H[lab]->GetMaximum());
    }

    const std::string cname = std::string("c_cnt_") + p.tag;
    TCanvas c(cname.c_str(), cname.c_str(), 900, 700);
    TH1D* frame = nullptr; for (int lab = 0; lab < nlabels; ++lab) if (H[lab]) { frame = H[lab].get(); break; }
    if (!frame) continue;
    frame->SetTitle((std::string(";") + (char)toupper(*p.tag) + "-plane semantic COUNT; " + (normalize_to_pdf ? "Probability density" : "Events")).c_str());
    frame->SetMaximum(1.3*ymax); frame->SetMinimum(0.0);
    frame->Draw("HIST");
    for (int lab = 0; lab < nlabels; ++lab) if (H[lab].get() != frame) H[lab]->Draw("HIST SAME");

    TLegend leg(0.60, 0.66, 0.88, 0.88);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42);
    for (int lab = 0; lab < nlabels; ++lab) if (H[lab]) leg.AddEntry(H[lab].get(), label_name(lab).c_str(), "l");
    leg.Draw();

    TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.04);
    tl.DrawLatex(0.14, 0.92, Form("#muBooNE Simulation – %s counts – Empty selection – %s",
                                  use_event_counts ? "event" : "slice",
                                  normalize_to_pdf ? "PDF" : "counts"));

    c.SaveAs((out_dir + "/semantic_counts_" + p.tag + ".pdf").c_str());
  }

  std::cout << "Saved all-label semantic COUNT PDFs (no channel split) in: " << out_dir << std::endl;
}
