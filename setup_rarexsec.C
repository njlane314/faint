#include "TInterpreter.h"
#include "TSystem.h"
#include "TError.h"
#include "TROOT.h"
#include "TString.h"
#include <cstring>
#include <string>

namespace {
static void load_one_macro(const TString& path) {
    if (path.IsNull()) return;
    if (gSystem->AccessPathName(path.Data())) return;
    const std::string inc = std::string("#include \"") + path.Data() + "\"";
    gInterpreter->ProcessLine(inc.c_str());
}
static void auto_load_macros(const char* incdir) {
    TString top = gSystem->Getenv("RAREXSEC");
    if (top.IsNull() && incdir && *incdir) {
        TString inc(incdir);
        if (inc.EndsWith("/rarexsec")) inc = gSystem->DirName(inc);
        if (inc.EndsWith("/include")) top = gSystem->DirName(inc);
        else top = gSystem->DirName(inc);
    }
    if (!top.IsNull()) {
        load_one_macro(top + "/scripts/rx_macros.C");
        load_one_macro(top + "/macros/rx_macros.C");
        load_one_macro(top + "/analysis/rx_macros.C");
    } else {
        load_one_macro("scripts/rx_macros.C");
        load_one_macro("macros/rx_macros.C");
        load_one_macro("analysis/rx_macros.C");
    }
}
}
namespace {
    static bool looks_like_macro(const std::string& s) {
        auto ends_with = [](const std::string& str, const char* suffix) {
            const std::size_t len = std::strlen(suffix);
            return str.size() >= len && str.compare(str.size() - len, len, suffix) == 0;
        };
        return ends_with(s, ".C") || ends_with(s, ".C+") || ends_with(s, ".C++") ||
               ends_with(s, ".cxx") || ends_with(s, ".cxx+") || ends_with(s, ".cxx++");
    }
}

void rx_call(const char* fname) {
    if (!fname || !*fname) { ::Error("rx_call","no function name"); return; }

    std::string code(fname);
    if (code.find('(') != std::string::npos) {
        gInterpreter->ProcessLine(code.c_str());
        return;
    }

    if (looks_like_macro(code)) {
        gROOT->Macro(code.c_str());
        return;
    }

    gInterpreter->ProcessLine((code + "()").c_str());
}
void setup_rarexsec(const char* libpath, const char* incdir) {
    if (libpath && *libpath) {
        const Int_t status = gSystem->Load(libpath);
        if (status < 0) { ::Error("setup_rarexsec","failed to load %s (%d)", libpath, status); }
    }
    if (incdir && *incdir) {
        gInterpreter->AddIncludePath(incdir);
        TString inc(incdir);
        if (inc.EndsWith("/rarexsec")) {
            TString parent = gSystem->DirName(inc);
            gInterpreter->AddIncludePath(parent.Data());
        } else {
            TString sub = inc + "/rarexsec";
            if (!gSystem->AccessPathName(sub.Data())) gInterpreter->AddIncludePath(sub.Data());
        }
    }
    auto_load_macros(incdir);
}
