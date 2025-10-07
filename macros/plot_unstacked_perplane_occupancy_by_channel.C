// plot_unstacked_perplane_occupancy_by_channel.C
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>
#include <TH1D.h>

#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>

#include <rarexsec/Hub.hh>
#include <rarexsec/Selection.hh>
#include <rarexsec/Plotter.hh>               // H1Spec, Options
#include <rarexsec/plot/UnstackedHist.hh>    // overlay plotter

// ---------- helpers ----------
static void load_libs(const char* extra_libs) {
  // Always try librexsec first; ignore if already loaded
  gSystem->Load("librarexsec");
  if (!extra_libs) return;
  std::string s(extra_libs);
  std::stringstream ss(s);
  std::string tok;
  auto try_load = [](const std::string& lib){
    if (lib.empty()) return;
    if (gSystem->Load(lib.c_str()) < 0)
      std::cerr << "Warning: failed to load '" << lib << "'\n";
  };
  // split by comma/semicolon/space
  for (char& c : s) if (c==';' || c==' ') c=',';
  std::stringstream ss2(s);
  while (std::getline(ss2, tok, ',')) try_load(tok);
}

// Fraction of pixels with ADC > thr in one image (vector<float>)
double occ_adc(const std::vector<float>& img, double thr) {
  if (img.empty()) return 0.0;
  std::size_t k = 0; for (float q : img) if (std::isfinite(q) && q > thr) ++k;
  return static_cast<double>(k) / static_cast<double>(img.size());
}

// ---------- main ----------
void plot_unstacked_perplane_occupancy_by_channel(const char* extra_libs = "",
                                                  double adc_thr = 0.0,
                                                  bool use_event_images = true,
                                                  bool normalize_to_pdf = true,
                                                  int nbins = 60,
                                                  double xmin = 0.0,
                                                  double xmax = 0.40,
                                                  const char* beamline = "numi-fhc",
                                                  const char* period   = "run1",
                                                  const char* config_path = "data/samples.json") {
  ROOT::EnableImplicitMT();
  load_libs(extra_libs);

  rarexsec::Hub hub(config_path);
  const std::vector<std::string> periods = {period};
  const auto mc = hub.simulation_entries(beamline, periods);
  if (mc.empty()) throw std::runtime_error("No MC samples found");

  const std::string out_dir = std::string("plots/") + beamline + "/" + period + "/unstacked_perplane_occ";
  gSystem->mkdir(out_dir.c_str(), true);

  // Common options
  rarexsec::plot::Options opt;
  opt.out_dir        = out_dir;
  opt.image_format   = "pdf";
  opt.show_ratio     = false;
  opt.legend_on_top  = true;
  opt.legend_split   = 0.85;
  opt.use_log_y      = false;
  opt.y_min          = 0.0;
  opt.y_max          = -1.0;
  opt.beamline       = beamline;
  opt.periods        = periods;

  // Choose image columns
  const char* Uimg = use_event_images ? "event_detector_image_u" : "detector_image_u";
  const char* Vimg = use_event_images ? "event_detector_image_v" : "detector_image_v";
  const char* Wimg = use_event_images ? "event_detector_image_w" : "detector_image_w";

  struct Plane { const char* tag; const char* col; };
  const Plane planes[3] = { {"u", Uimg}, {"v", Vimg}, {"w", Wimg} };

  for (const auto& p : planes) {
    rarexsec::plot::H1Spec spec{
      .id     = std::string("image_occupancy_") + p.tag,
      .title  = (std::string(";") + (char)toupper(*p.tag) + "-plane image occupancy;Events"),
      .expr   = std::string("occ_adc(") + p.col + "," + std::to_string(adc_thr) + ")",
      .weight = "w_nominal",
      .nbins  = nbins,
      .xmin   = xmin,
      .xmax   = xmax,
      .sel    = rarexsec::selection::Preset::Empty
    };

    opt.x_title = std::string(1, (char)toupper(*p.tag)) + "-plane image occupancy";
    opt.y_title = "Events";

    rarexsec::plot::UnstackedHist plot(spec, opt, mc, /*data*/{}, normalize_to_pdf, /*line_width*/3);
    plot.draw_and_save("pdf");
  }

  std::cout << "Saved per-plane unstacked occupancy PDFs to: " << out_dir << std::endl;
}
