// scripts/setup_faint.C
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "TError.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TSystem.h"

namespace {

class ErrorLevelGuard {
public:
  explicit ErrorLevelGuard(int new_level)
      : previous_level_(gErrorIgnoreLevel) {
    gErrorIgnoreLevel = new_level;
  }

  ErrorLevelGuard(const ErrorLevelGuard&) = delete;
  ErrorLevelGuard& operator=(const ErrorLevelGuard&) = delete;

  ~ErrorLevelGuard() { gErrorIgnoreLevel = previous_level_; }

private:
  int previous_level_;
};

std::string get_env(const char* name) {
  if (name == nullptr) {
    return {};
  }
  if (const char* value = gSystem->Getenv(name)) {
    return value;
  }
  return {};
}

std::string parent_dir(std::string path) {
  if (path.empty()) {
    return {};
  }

  // Strip trailing separators to avoid returning the same directory.
  while (!path.empty() && (path.back() == '/' || path.back() == '\\')) {
    path.pop_back();
  }

  const std::string::size_type pos = path.find_last_of("/\\");
  if (pos == std::string::npos) {
    return {};
  }
  if (pos == 0) {
    return path.substr(0, 1);
  }
  return path.substr(0, pos);
}

std::string infer_topdir() {
  // Allow users to point to a custom installation prefix.
  std::string from_env = get_env("FAINT");
  if (!from_env.empty()) {
    return from_env;
  }
  from_env = get_env("FAINT_PREFIX");
  if (!from_env.empty()) {
    return from_env;
  }

  // `Which` searches the ROOT macro path for this helper and returns the
  // absolute location if it is discoverable.
  if (char* located = gSystem->Which(gROOT->GetMacroPath(), "setup_faint.C")) {
    std::string macro_path = located;
    std::free(located);

    const std::string macro_dir = parent_dir(macro_path);
    if (!macro_dir.empty()) {
      return parent_dir(macro_dir);
    }
  }

  return {};
}

std::string join_path(const std::string& base, const std::string& leaf) {
  if (base.empty()) {
    return {};
  }
  if (!leaf.empty() && (leaf.front() == '/' || leaf.front() == '\\')) {
    return leaf;
  }

  if (base.back() == '/' || base.back() == '\\') {
    return base + leaf;
  }
  return base + "/" + leaf;
}

} // namespace
void load_header(const std::string& h) {
  TInterpreter::EErrorCode ec;
  gInterpreter->ProcessLine( (std::string("#include <") + h + ">").c_str(), &ec );
  if (ec != 0) {
    std::cout << "Error including header <" << h << ">. "
              << "Ensure your include directory is in ROOT_INCLUDE_PATH.\n";
  }
}

void setup_faint(const char* abs_lib_path = nullptr, const char* abs_inc_dir = nullptr) {
  ErrorLevelGuard error_level_guard(kFatal);
  // Some ROOT 6 builds need libGraf preloaded for dictionaries
  if (gROOT->GetVersionInt() >= 60000) {
    if (gSystem->Load("libGraf") != 0) {
      std::cout << "Warning: could not preload libGraf\n";
    }
  }

  const std::string topdir = infer_topdir();

  std::string include_dir;
  if (abs_inc_dir && std::strlen(abs_inc_dir) > 0) {
    include_dir = abs_inc_dir;
  } else {
    include_dir = join_path(topdir, "include");
  }
  if (!include_dir.empty()) {
    gSystem->AddIncludePath( (std::string("-I") + include_dir).c_str() );
  } else {
    std::cout << "Warning: could not determine the faint include directory. "
              << "Pass it explicitly to setup_faint().\n";
  }

  // Load library (absolute path if given, else rely on DYLD/LD_LIBRARY_PATH)
  int rc = 1;
  std::string lib_path;
  bool attempted_absolute = false;
  if (abs_lib_path && std::strlen(abs_lib_path) > 0) {
    lib_path = abs_lib_path;
    attempted_absolute = true;
    rc = gSystem->Load(lib_path.c_str());
  } else {
#ifdef __APPLE__
    lib_path = join_path(topdir, "build/lib/libfaint.dylib");
#else
    lib_path = join_path(topdir, "build/lib/libfaint.so");
#endif
    if (!lib_path.empty()) {
      attempted_absolute = true;
      rc = gSystem->Load(lib_path.c_str());
    }
  }
  if (rc != 0 && attempted_absolute && !lib_path.empty()) {
    std::cout << "Warning: failed to load faint library from '"
              << lib_path
              << "'.\n";
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
  if (!include_dir.empty()) {
    load_header("faint/Types.h");
    load_header("faint/SampleSet.h");
    load_header("faint/PreSelection.h");
    load_header("faint/TruthClassifier.h");
    load_header("faint/MuonSelector.h");
    load_header("faint/Weighter.h");
  }
}
