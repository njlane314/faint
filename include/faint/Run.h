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
  RunReader() = default;
  explicit RunReader(const std::string& path);

  void add(Run rc);

  const Run& get(const std::string& beam, const std::string& period) const;

  const std::map<std::string, Run>& all() const noexcept;

private:
  void load_from_json(const nlohmann::json& data);

  std::map<std::string, Run> configs_;
};

}  // namespace faint

#endif
