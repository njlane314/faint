#pragma once
#include <ROOT/RDataFrame.hxx>
#include <string>

namespace rarexsec {

class Processor {
public:
    ROOT::RDF::RNode run(ROOT::RDF::RNode node,
                         sample::origin kind) const;
};

const Processor& processor();

} 