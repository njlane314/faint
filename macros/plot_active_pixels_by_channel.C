// plot_active_pixels_by_channel.C
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>
#include <TSystem.h>
#include <TH1D.h>

#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>

#include <rarexsec/Hub.hh>
#include <rarexsec/proc/Selection.hh>
#include <rarexsec/Plotter.hh>               // H1Spec, Options, Plotter

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

// Return the ABSOLUTE COUNT of active pixels (> thr) in a single image container
namespace {

template <typename Image>
double active_pixels_impl(const Image& img, double thr) {
  if (img.empty()) return 0.0;
  std::size_t k = 0; for (float q : img) if (std::isfinite(q) && q > thr) ++k;
  return static_cast<double>(k);
}

}  // namespace

double active_pixels(const std::vector<float>& img, double thr) {
  return active_pixels_impl(img, thr);
}

// Overload for RDF columns backed by RVec<float>
double active_pixels(const ROOT::VecOps::RVec<float>& img, double thr) {
  return active_pixels_impl(img, thr);
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

  const bool has_xmax_override = xmax > 0.0;
  const double xmax_override   = xmax;
  const double NPIX_plane      = 512.0 * 512.0;            // per plane

  struct PlotConfig {
    std::string id;
    std::string title;
    std::string expr;
    std::string x_title;
    double default_xmax;
  };

  std::vector<PlotConfig> plot_configs;
  plot_configs.reserve(3);

  const auto make_expr = [adc_thr](const char* column) {
    return std::string("active_pixels(") + column + "," + std::to_string(adc_thr) + ")";
  };

  if (plane_upper == "U") {
    plot_configs.push_back({
      .id            = "active_pixels_u",
      .title         = ";Active pixels (U plane);Events",
      .expr          = make_expr(Uimg),
      .x_title       = "Active pixels (U plane)",
      .default_xmax  = NPIX_plane
    });
  } else if (plane_upper == "V") {
    plot_configs.push_back({
      .id            = "active_pixels_v",
      .title         = ";Active pixels (V plane);Events",
      .expr          = make_expr(Vimg),
      .x_title       = "Active pixels (V plane)",
      .default_xmax  = NPIX_plane
    });
  } else if (plane_upper == "W") {
    plot_configs.push_back({
      .id            = "active_pixels_w",
      .title         = ";Active pixels (W plane);Events",
      .expr          = make_expr(Wimg),
      .x_title       = "Active pixels (W plane)",
      .default_xmax  = NPIX_plane
    });
  } else if (plane_upper == "SUM" || plane_upper == "ALL" || plane_upper == "UVW" || plane_upper.empty()) {
    plot_configs.push_back({
      .id            = "active_pixels_u",
      .title         = ";Active pixels (U plane);Events",
      .expr          = make_expr(Uimg),
      .x_title       = "Active pixels (U plane)",
      .default_xmax  = NPIX_plane
    });
    plot_configs.push_back({
      .id            = "active_pixels_v",
      .title         = ";Active pixels (V plane);Events",
      .expr          = make_expr(Vimg),
      .x_title       = "Active pixels (V plane)",
      .default_xmax  = NPIX_plane
    });
    plot_configs.push_back({
      .id            = "active_pixels_w",
      .title         = ";Active pixels (W plane);Events",
      .expr          = make_expr(Wimg),
      .x_title       = "Active pixels (W plane)",
      .default_xmax  = NPIX_plane
    });
  } else {
    throw std::runtime_error("Unknown plane: " + plane_key);
  }

  rarexsec::plot::Options base_opt;
  base_opt.out_dir        = out_dir;
  base_opt.image_format   = "pdf";
  base_opt.show_ratio     = false;            // no ratio panel for raw counts here
  base_opt.legend_on_top  = true;
  base_opt.legend_split   = 0.85;
  base_opt.use_log_y      = false;            // turn on if small tails are hard to see
  base_opt.y_min          = 0.0;
  base_opt.y_max          = -1.0;
  base_opt.beamline       = beamline;
  base_opt.periods        = periods;
  base_opt.y_title        = "Events";

  for (const auto& cfg : plot_configs) {
    rarexsec::plot::H1Spec spec{
      .id     = cfg.id,
      .title  = cfg.title,
      .expr   = cfg.expr,              // evaluated via RDF Define in UnstackedHist
      .weight = "w_nominal",
      .nbins  = nbins,
      .xmin   = xmin,
      .xmax   = has_xmax_override ? xmax_override : cfg.default_xmax,
      .sel    = rarexsec::selection::Preset::Empty
    };

    auto opt = base_opt;
    opt.x_title = cfg.x_title;

    // Draw overlay by analysis channel (NO normalization)
    rarexsec::plot::Plotter plotter(opt);
    plotter.draw_unstacked_by_channel(spec, mc, /*normalize_to_pdf*/ false, /*line_width*/3);

    std::cout << "Saved: " << out_dir << "/" << cfg.id << ".pdf" << std::endl;
  }
}
