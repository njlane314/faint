#ifndef FAINT_RUN_CATALOG_H
#define FAINT_RUN_CATALOG_H

#include <map>
#include <string>

#include "nlohmann/json.hpp"

#include "faint/Run.h"

namespace faint {

class RunCatalog {
 public:
  RunCatalog() = default;

  const Run& get(const std::string& beam, const std::string& period) const;

  const std::map<std::string, Run>& all() const noexcept { return runs_; }

  bool empty() const noexcept { return runs_.empty(); }

  static RunCatalog from_json(const nlohmann::json& data);

  static RunCatalog from_file(const std::string& path);

 private:
  std::map<std::string, Run> runs_;
};

}  // namespace faint

#endif  // FAINT_RUN_CATALOG_H

