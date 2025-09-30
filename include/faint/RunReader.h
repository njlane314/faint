#ifndef FAINT_RUN_READER_H
#define FAINT_RUN_READER_H

#include <map>
#include <string>

#include "nlohmann/json.hpp"

#include <faint/Run.h>

namespace faint {

class RunReader {
public:
  void add(Run rc);

  const Run& get(const std::string& beam, const std::string& period) const;

  const std::map<std::string, Run>& all() const noexcept;

  static RunReader from_json(const nlohmann::json& data);

  static RunReader from_file(const std::string& path);

  static RunReader read(const nlohmann::json& data) { return from_json(data); }

  static RunReader read(const std::string& path) { return from_file(path); }

private:
  std::map<std::string, Run> configs_;
};

}

#endif
