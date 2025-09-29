#ifndef RUN_CONFIG_H
#define RUN_CONFIG_H

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

}

#endif
