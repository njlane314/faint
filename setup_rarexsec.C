#include "TInterpreter.h"
#include "TSystem.h"
#include "TError.h"

#include <string>

namespace {

const char* const kHeadersToInclude[] = {
    "rarexsec/Hub.hh",
    "rarexsec/Processor.hh",
    "rarexsec/Plotter.hh",
    "rarexsec/proc/Selection.hh",
    "rarexsec/proc/Snapshot.hh",
    "rarexsec/proc/Volume.hh",
    "rarexsec/plot/StackedHist.hh",
    "rarexsec/plot/Channels.hh"
};

void include_required_headers() {
    for (const auto* header : kHeadersToInclude) {
        const std::string directive = std::string("#include \"") + header + "\"";
        gInterpreter->ProcessLine(directive.c_str());
    }
}

} 

void setup_rarexsec(const char* libpath, const char* incdir) {
    if (libpath && libpath[0] != '\0') {
        const Int_t loadStatus = gSystem->Load(libpath);
        if (loadStatus < 0) {
            ::Error("setup_rarexsec", "Failed to load library '%s' (status %d)", libpath, loadStatus);
        }
    }

    if (incdir && incdir[0] != '\0') {
        gInterpreter->AddIncludePath(incdir);
    }

    include_required_headers();
}
