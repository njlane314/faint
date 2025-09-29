#ifndef SAMPLE_DATASET_H
#define SAMPLE_DATASET_H

#include <unordered_map>

#include "ROOT/RDataFrame.hxx"

#include <rarexsec/core/AnalysisKey.h>
#include <rarexsec/data/SampleTypes.h>

namespace analysis {

struct Dataset {
    SampleOrigin origin_;
    AnalysisRole role_;
    mutable ROOT::RDF::RNode dataframe_;
};

struct DatasetVariations {
    Dataset nominal_;
    std::unordered_map<SampleVariation, Dataset> variations_;
};

using DatasetMap = std::unordered_map<SampleKey, DatasetVariations>;

}

#endif
