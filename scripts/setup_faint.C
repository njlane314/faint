// scripts/setup_faint.C
#include <cstring>
#include <iostream>
#include <string>

#include "TInterpreter.h"
#include "TROOT.h"
#include "TSystem.h"
void load_header(const std::string& h) {
  TInterpreter::EErrorCode ec;
  gInterpreter->ProcessLine( (std::string("#include <") + h + ">").c_str(), &ec );
  if (ec != 0) {
    std::cout << "Error including header <" << h << ">. "
              << "Ensure your include directory is in ROOT_INCLUDE_PATH.\n";
  }
}

void setup_faint(const char* abs_lib_path = nullptr, const char* abs_inc_dir = nullptr) {
  // Some ROOT 6 builds need libGraf preloaded for dictionaries
  if (gROOT->GetVersionInt() >= 60000) {
    if (gSystem->Load("libGraf") != 0) {
      std::cout << "Warning: could not preload libGraf\n";
    }
  }

  // Include dir (optional explicit arg)
  if (abs_inc_dir && std::strlen(abs_inc_dir) > 0) {
    gSystem->AddIncludePath( (std::string("-I") + abs_inc_dir).c_str() );
  }

  // Load library (absolute path if given, else rely on DYLD/LD_LIBRARY_PATH)
  int rc = 1;
  if (abs_lib_path && std::strlen(abs_lib_path) > 0) {
    rc = gSystem->Load(abs_lib_path);
  }
  if (rc != 0) {
    rc = gSystem->Load("libfaint"); // fall back to soname
  }

  if (rc == 0) {
    std::cout << "Loaded faint library.\n";
  } else {
    const char* env_var = gSystem->UnixPathName("/") ? "LD_LIBRARY_PATH" : "DYLD_LIBRARY_PATH";
    std::cout << "Error loading faint library. Add its directory to "
              << env_var
              << " or pass an absolute path to setup_faint().\n";
  }

  // Pull in commonly used headers so macros can use the API immediately
  load_header("faint/Types.h");
  load_header("faint/SampleSet.h");
  load_header("faint/NuMuCCSelector.h");
  load_header("faint/TruthClassifier.h");
  load_header("faint/MuonSelector.h");
  load_header("faint/Weighter.h");
}
