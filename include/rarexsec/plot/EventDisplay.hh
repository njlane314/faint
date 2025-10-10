#ifndef RAREXSEC_PLOT_EVENTDISPLAY_HH
#define RAREXSEC_PLOT_EVENTDISPLAY_HH

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"
#include "TStyle.h"

#include <ROOT/RDataFrame.hxx>

namespace rarexsec {
namespace plot {

class EventDisplay {
public:
    enum class Mode { Detector, Semantic };

    static Mode parse_mode(const std::string& s) {
        if (s == "semantic" || s == "Semantic") return Mode::Semantic;
        return Mode::Detector;
    }

    struct Spec {
        std::string id;          ///< Unique identifier (used in file names)
        std::string title;       ///< Canvas / histogram title
        Mode        mode{Mode::Detector};
        int         grid_w{0};   ///< Optional; if 0, deduced assuming square
        int         grid_h{0};   ///< Optional; if 0, deduced assuming square
    };

    struct Options {
        std::string out_dir = "plots";
        int         canvas_size = 800;
        double      margin = 0.10;      ///< symmetric margins (fraction)
        bool        use_log_z = true;   ///< detector mode only

        // ---- Detector mode tweaks
        double      det_threshold = 4.0;
        double      det_min       = 1.0;
        double      det_max       = 1000.0;

        // ---- Semantic mode tweaks
        bool        show_legend   = true;
        int         legend_cols   = 5;  ///< semantic legend columns
    };

    using DetectorData = std::vector<float>;
    using SemanticData = std::vector<int>;

    // Constructors for detector / semantic images
    EventDisplay(Spec spec, Options opt, DetectorData data);
    EventDisplay(Spec spec, Options opt, SemanticData data);

    // Draw this one display
    void draw(TCanvas& canvas);

    // Draw and save (auto-creates canvas like StackedHist::draw_and_save)
    void draw_and_save(const std::string& image_format = "png");
    void draw_and_save(const std::string& image_format, const std::string& file_override);

    // ----------------- OPTIONAL: batched rendering from RDataFrame -----------------
    // Still only ONE class (this one). The batch options live inside the class.
    struct BatchOptions {
        // ---- DataFrame selection
        std::string selection_expr;       // empty => no filter
        unsigned long long n_events{1};

        // ---- I/O
        std::string out_dir{"./plots/event_displays"};
        std::string image_format{"png"};
        std::string combined_pdf;         // if non-empty and image_format == "pdf"
        std::string manifest_path;        // if non-empty, write a json manifest

        // ---- Planes & columns
        std::vector<std::string> planes{"U","V","W"};

        struct Columns {
            std::string run = "run";
            std::string sub = "sub";
            std::string evt = "evt";
            std::string det_u = "event_detector_image_u";
            std::string det_v = "event_detector_image_v";
            std::string det_w = "event_detector_image_w";
            std::string sem_u = "semantic_image_u";
            std::string sem_v = "semantic_image_v";
            std::string sem_w = "semantic_image_w";
        } cols;

        // ---- File naming
        std::string file_pattern{"{plane}_{run}_{sub}_{evt}"};

        // ---- Mode & per-image display options
        Mode     mode{Mode::Detector};
        Options  display; // canvas size, margins, logz, legend, etc.
    };

    // Render event displays from an RDF node according to BatchOptions.
    // If 'combined_pdf' is set and 'image_format'=="pdf", output is a single
    // multi-page PDF using ROOT's "file(", "file", "file)" protocol.
    static void render_from_rdf(ROOT::RDF::RNode df, const BatchOptions& opt);

private:
    // lifecycle helpers
    void setup_canvas(TCanvas& c) const;
    void build_histogram();  // creates TH2 from data based on mode

    // draw paths
    void draw_detector(TCanvas& c);
    void draw_semantic(TCanvas& c);
    void draw_semantic_legend();

    // utilities
    static std::pair<int,int> deduce_grid(int requested_w,
                                          int requested_h,
                                          std::size_t flat_size);

private:
    Spec  spec_;
    Options opt_;

    std::variant<DetectorData, SemanticData> data_;

    std::unique_ptr<TH2F>    hist_;
    std::unique_ptr<TLegend> legend_;
    std::vector<std::unique_ptr<TH1F>> legend_entries_;

    std::string plot_name_;
    std::string output_directory_;
};

} // namespace plot
} // namespace rarexsec

#endif // RAREXSEC_PLOT_EVENTDISPLAY_HH
