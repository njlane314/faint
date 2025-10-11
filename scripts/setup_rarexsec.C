#include "TSystem.h"
#include "TInterpreter.h"
#include <cstdio>
#include <string>

void setup_rarexsec(const char* libpath, const char* incdir) {
    if (libpath && *libpath) {
        int s = gSystem->Load(libpath);
        if (s < 0) std::fprintf(stderr, "setup_rarexsec: failed to load %s (%d)\n", libpath, s);
    }

    if (incdir && *incdir) {
        gInterpreter->AddIncludePath(incdir);
        std::string sub = std::string(incdir) + "/rarexsec";
        if (!gSystem->AccessPathName(sub.c_str())) gInterpreter->AddIncludePath(sub.c_str());
    }
}
