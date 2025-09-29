#ifndef FAINT_RUN_CATALOG_H
#define FAINT_RUN_CATALOG_H

#include <map>
#include <string>
#include <utility>

#include "nlohmann/json.hpp"

#include "faint/Run.h"
#include "faint/RunReader.h"

namespace faint {

class RunCatalog {
 public:
  RunCatalog() = default;

  explicit RunCatalog(RunReader reader) : reader_(std::move(reader)) {}

  const Run& get(const std::string& beam, const std::string& period) const {
    return reader_.get(beam, period);
  }

  const std::map<std::string, Run>& all() const noexcept { return reader_.all(); }

  static RunCatalog fromJson(const nlohmann::json& data) {
    return RunCatalog{RunReader::read(data)};
  }

  static RunCatalog fromFile(const std::string& path) {
    return RunCatalog{RunReader::read(path)};
  }

 private:
  RunReader reader_;
};

}  // namespace faint

#endif
