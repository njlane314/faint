#pragma once
#include <ROOT/RDataFrame.hxx>
#include <string>

namespace rarexsec {

struct Entry;

class Processor {
public:
    ROOT::RDF::RNode run(ROOT::RDF::RNode node,
                         const Entry& rec) const;
};

const Processor& processor();

} 