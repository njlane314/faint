#ifndef FAINT_CAMPAIGN_H
#define FAINT_CAMPAIGN_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <faint/study.h>

namespace analysis {
namespace study {

class Campaign {
public:
    Campaign() = default;

    explicit Campaign(Study study)
        : study_(std::make_shared<Study>(std::move(study))) {}

    static Campaign open(const std::string& run_config_json, Options opt, Variables vars = Variables{}) {
        return Campaign{Study::open(run_config_json, std::move(opt), std::move(vars))};
    }

    bool has_study() const noexcept { return static_cast<bool>(study_); }

    const Study& study() const {
        if (!study_) throw std::runtime_error("Campaign study has not been initialised");
        return *study_;
    }

    Study& study() {
        if (!study_) throw std::runtime_error("Campaign study has not been initialised");
        return *study_;
    }

    void set_study(Study study) {
        study_ = std::make_shared<Study>(std::move(study));
    }

private:
    std::shared_ptr<Study> study_;
};

} // namespace study
} // namespace analysis

#endif
