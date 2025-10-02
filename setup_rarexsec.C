#include "TInterpreter.h"
#include "TSystem.h"
#include "TError.h"

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
}
