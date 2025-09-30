#include "faint/Selection.h"

#include <utility>

namespace faint {

Selection::Selection() = default;

Selection::Selection(std::string expression)
    : expression_(std::move(expression)) {}

const std::string& Selection::str() const noexcept { return expression_; }

bool Selection::empty() const noexcept { return expression_.empty(); }

}  // namespace faint
