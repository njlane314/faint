#pragma once
#include <string>
#include <string_view>
#include <vector>

namespace rarexsec {

struct Entry;
class Hub;

namespace snapshot {

struct Options {
    std::string outdir = "snapshots";
    std::string tree = "analysis";
    std::vector<std::string> columns;
};

std::vector<std::string> write(const std::vector<const Entry*>& samples,
                               const Options& opt = {});

std::vector<std::string> write(const Hub& hub,
                               std::string_view beamline,
                               const std::vector<std::string>& periods,
                               const Options& opt = {});

}
}
