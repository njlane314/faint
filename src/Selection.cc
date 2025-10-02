#include "rarexsec/Selection.hh"

namespace rarexsec {
namespace selection {

bool passes_muon_track_selection(float score,
                                 float llr,
                                 float length,
                                 float distance,
                                 unsigned generation) {
    return score > min_score &&
           llr > min_llr &&
           length > min_length &&
           distance < max_distance &&
           generation == required_generation;
}

}  // namespace selection
}  // namespace rarexsec
