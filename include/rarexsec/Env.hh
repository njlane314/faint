#pragma once

#include <string>
#include <vector>

#include "rarexsec/Hub.hh"
#include "rarexsec/Processor.hh"

namespace rarexsec {
struct Env {
  std::string cfg, tree, beamline;
  std::vector<std::string> periods;
  static Env from_env();
  Hub make_hub() const {
    ProcessorOptions o;
    o.tree = tree;
    return Hub(cfg, o);
  }
};
} // namespace rarexsec
