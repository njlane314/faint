#include <ROOT/RDataFrame.hxx>
#include <faint/Dataset.h>
#include <faint/Log.h>

#include <iostream>

int main(int argc, char** argv) {
  faint::log::init();
  faint::log::set_level(faint::log::Level::kDebug);

  std::string config_path;
  if (argc > 1) {
    config_path = argv[1];
  } else {
    config_path = faint::dataset::run_config_path();
  }

  faint::dataset::Options options;
  options.ntuple_dir = faint::dataset::ntuple_directory(config_path);
  auto dataset = faint::dataset::Dataset::open(config_path, options);

  auto keys = dataset.sample_keys();
  std::cout << "Loaded the following samples:\n";
  for (const auto& key : keys) {
    std::cout << "  - " << key << '\n';
  }

  auto df = dataset.df(keys.at(0));
  auto count = df.Count();
  std::cout << "Count of rows: " << *count << '\n';

  return 0;
}
