#ifndef FAINT_RUN_H
#define FAINT_RUN_H

#include <map>
#include <string>

#include <nlohmann/json.hpp>

namespace faint {

using json = nlohmann::json;

struct Run {
    std::string beam_mode;
    std::string run_period;

    double nominal_pot{0.0};
    long nominal_triggers{0};

    json samples;

    Run(const json& j, std::string bm, std::string pr);

    std::string key() const;
    std::string label() const;

    void validate() const;
};

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

}  // namespace faint

#endif
