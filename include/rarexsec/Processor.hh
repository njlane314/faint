#pragma once
#include "rarexsec/proc/DataModel.hh"
#include <ROOT/RDataFrame.hxx>

namespace rarexsec {

class Processor {
  public:
    ROOT::RDF::RNode run(ROOT::RDF::RNode node, const rarexsec::Entry& rec) const;
};

const Processor& processor();

}
