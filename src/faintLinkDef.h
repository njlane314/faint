#include "faint/Dataset.h"
#include "faint/proc/EventProcessor.h"
#include "faint/proc/MuonSelector.h"
#include "faint/proc/PreSelection.h"
#include "faint/Run.h"
#include "faint/Samples.h"
#include "faint/Selection.h"
#include "faint/plot/StackedHistogram.h"
#include "faint/proc/TruthClassifier.h"
#include "faint/Variables.h"
#include "faint/proc/Weighter.h"

#ifdef __CLING__
#pragma link C++ namespace faint;
#pragma link C++ namespace faint::dataset;
#pragma link C++ namespace faint::plot;

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;

#pragma link C++ class faint::sample::Key+;
#pragma link C++ enum faint::sample::Origin;
#pragma link C++ enum faint::sample::Role;
#pragma link C++ enum faint::sample::Variation;

#pragma link C++ class faint::dataset::Options+;
#pragma link C++ class faint::dataset::Dataset+;
#pragma link C++ class faint::dataset::Dataset::Entry+;
#pragma link C++ class faint::dataset::Dataset::Variations+;

#pragma link C++ class faint::Sample+;
#pragma link C++ class faint::SampleSet+;
#pragma link C++ class faint::Run+;
#pragma link C++ class faint::RunReader+;
#pragma link C++ class faint::Selection+;
#pragma link C++ class faint::Variables+;
#pragma link C++ class faint::Weighter+;
#pragma link C++ class faint::plot::StackedHistogram+;

#pragma link C++ class faint::EventProcessor+;
#pragma link C++ class faint::PreSelection+;
#pragma link C++ class faint::MuonSelector+;
#pragma link C++ class faint::TruthClassifier+;
#endif
