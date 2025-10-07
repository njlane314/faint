// plot_active_pixels_by_channel.C
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>
#include <TH1D.h>

#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

#include <rarexsec/Hub.hh>
#include <rarexsec/proc/Selection.hh>
#include <rarexsec/Plotter.hh>               // H1Spec, Options
#include <rarexsec/plot/UnstackedHist.hh>    // overlay-by-channel plotter

// ---------------- Helpers ----------------

// Load librexsec and any additional user libs (comma/space/semicolon separated)
static void load_libs(const char* extra_libs) {
  gSystem->Load("librarexsec");
  if (!extra_libs) return;
  std::string s(extra_libs);
  for (char& c : s) if (c==';' || c==' ') c=',';
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) if (!tok.empty()) {
    if (gSystem->Load(tok.c_str()) < 0)
      std::cerr << "Warning: failed to load '" << tok << "'\n";
  }
}

// Return the ABSOLUTE COUNT of active pixels (> thr) in a single image (vector<float>)
double active_pixels(const std::vector<float>& img, double thr) {
  if (img.empty()) return 0.0;
  std::size_t k = 0; for (float q : img) if (std::isfinite(q) && q > thr) ++k;
  return static_cast<double>(k);
}

// ---------------- Main entry ----------------
void plot_active_pixels_by_channel(const char* extra_libs   = "",
                                   const char* plane        = "sum",   // "u", "v", "w", or "sum"
                                   double adc_thr           = 4.0,     // default to your image threshold
                                   bool use_event_images    = true,    // false -> slice images
                                   int nbins                = 120,
                                   double xmin              = 0.0,
                                   double xmax              = -1.0,    // auto: NPIX (plane) or 3*NPIX (sum)
                                   const char* beamline     = "numi-fhc",
                                   const char* period       = "run1",
                                   const char* config_path  = "data/samples.json")
{
  ROOT::EnableImplicitMT();
  load_libs(extra_libs);

  // Inputs
  rarexsec::Hub hub(config_path);
  const std::vector<std::string> periods = {period};
  const auto mc = hub.simulation_entries(beamline, periods);
  if (mc.empty()) throw std::runtime_error("No MC samples found");

  // Output dir
  const std::string out_dir = std::string("plots/") + beamline + "/" + period + "/active_pixels";
  gSystem->mkdir(out_dir.c_str(), true);

  // Pick image columns
  const char* Uimg = use_event_images ? "event_detector_image_u" : "detector_image_u";
  const char* Vimg = use_event_images ? "event_detector_image_v" : "detector_image_v";
  const char* Wimg = use_event_images ? "event_detector_image_w" : "detector_image_w";

  // Axis range (absolute counts)
  auto up = [](char c){ return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); };
  const std::string plane_key = plane ? plane : "sum";
  std::string plane_upper;
  plane_upper.reserve(plane_key.size());
  for (char c : plane_key) plane_upper.push_back(up(c));

  const double NPIX_plane = 512.0 * 512.0;            // per plane
  const double NPIX_sum   = 3.0 * NPIX_plane;         // U+V+W
  if (xmax <= 0.0) {
    xmax = (plane_upper == "SUM") ? NPIX_sum : NPIX_plane;
  }

  // Build expression for the chosen plane/sum
  std::string expr, x_title, id, title;

  if (plane_upper == "U") {
    expr   = std::string("active_pixels(") + Uimg + "," + std::to_string(adc_thr) + ")";
    x_title= "Active pixels (U plane)";
    id     = "active_pixels_u";
    title  = ";Active pixels (U plane);Events";
  } else if (plane_upper == "V") {
    expr   = std::string("active_pixels(") + Vimg + "," + std::to_string(adc_thr) + ")";
    x_title= "Active pixels (V plane)";
    id     = "active_pixels_v";
    title  = ";Active pixels (V plane);Events";
  } else if (plane_upper == "W") {
    expr   = std::string("active_pixels(") + Wimg + "," + std::to_string(adc_thr) + ")";
    x_title= "Active pixels (W plane)";
    id     = "active_pixels_w";
    title  = ";Active pixels (W plane);Events";
  } else {
    // Sum across U+V+W
    expr   = std::string("active_pixels(") + Uimg + "," + std::to_string(adc_thr) + ")"
           + " + active_pixels(" + Vimg + "," + std::to_string(adc_thr) + ")"
           + " + active_pixels(" + Wimg + "," + std::to_string(adc_thr) + ")";
    x_title= "Active pixels (U+V+W sum)";
    id     = "active_pixels_sum";
    title  = ";Active pixels (U+V+W);Events";
  }

  // Plot spec and options
  rarexsec::plot::H1Spec spec{
    .id     = id,
    .title  = title,
    .expr   = expr,              // evaluated via RDF Define in UnstackedHist
    .weight = "w_nominal",
    .nbins  = nbins,
    .xmin   = xmin,
    .xmax   = xmax,
    .sel    = rarexsec::selection::Preset::Empty
  };

  rarexsec::plot::Options opt;
  opt.out_dir        = out_dir;
  opt.image_format   = "pdf";
  opt.show_ratio     = false;            // no ratio panel for raw counts here
  opt.legend_on_top  = true;
  opt.legend_split   = 0.85;
  opt.use_log_y      = false;            // turn on if small tails are hard to see
  opt.y_min          = 0.0;
  opt.y_max          = -1.0;
  opt.beamline       = beamline;
  opt.periods        = periods;
  opt.x_title        = x_title;
  opt.y_title        = "Events";

  // Draw overlay by analysis channel (NO normalization)
  rarexsec::plot::UnstackedHist plot(spec, opt, mc, /*data*/{}, /*normalize_to_pdf*/ false, /*line_width*/3);
  plot.draw_and_save("pdf");

  std::cout << "Saved: " << out_dir << "/" << id << ".pdf" << std::endl;
}
