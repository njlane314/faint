#pragma once
#include <ROOT/RDataFrame.hxx>
#include "rarexsec/DataModel.hh"

namespace rarexsec {

class Processor {
public:
  ROOT::RDF::RNode run(ROOT::RDF::RNode node,
                       const rarexsec::Entry& rec,
                       const ProcessorOptions& opt) const;
};

const Processor& processor();

} // namespace rarexsec
