#include <cstring>
#include <iostream>
#include <string>

#include "TInterpreter.h"
#include "TROOT.h"
#include "TSystem.h"

void load_header(const std::string& h) {
  TInterpreter::EErrorCode ec;
  gInterpreter->ProcessLine((std::string("#include <") + h + ">").c_str(), &ec);
  if (ec != 0) {
    std::cout << "Error including header <" << h << ">. "
              << "Ensure your include directory is in ROOT_INCLUDE_PATH.\n";
  }
}

void setup_rarexsec(const char* abs_lib_path = nullptr, const char* abs_inc_dir = nullptr) {
  // Some ROOT 6 builds need libGraf preloaded for dictionaries
  if (gROOT->GetVersionInt() >= 60000) {
    if (gSystem->Load("libGraf") != 0) {
      std::cout << "Warning: could not preload libGraf\n";
    }
  }

  // Include dir (optional explicit arg)
  if (abs_inc_dir && std::strlen(abs_inc_dir) > 0) {
    gSystem->AddIncludePath((std::string("-I") + abs_inc_dir).c_str());
  }

  // Load library (absolute path if given, else rely on DYLD/LD_LIBRARY_PATH)
  int rc = 1;
  if (abs_lib_path && std::strlen(abs_lib_path) > 0) {
    rc = gSystem->Load(abs_lib_path);
  }
  if (rc != 0) {
    rc = gSystem->Load("librarexsec"); // fall back to soname
  }

  if (rc == 0) {
    std::cout << "Loaded rarexsec library.\n";
  } else {
    std::cout << "Error loading rarexsec library. Add its directory to "
              << (gSystem->UnixPathName("/") ? "LD_LIBRARY_PATH" : "DYLD_LIBRARY_PATH")
              << " or pass an absolute path to setup_rarexsec().\n";
  }

  // Pull in commonly used headers so macros can use the API immediately
  load_header("rarexsec/data/Types.h");
  load_header("rarexsec/data/SampleSet.h");
  load_header("rarexsec/data/NuMuCCSelector.h");
  load_header("rarexsec/data/TruthClassifier.h");
  load_header("rarexsec/data/MuonSelector.h");
  load_header("rarexsec/data/Weighter.h");
}
