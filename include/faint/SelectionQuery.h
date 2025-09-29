#ifndef FAINT_SELECTION_QUERY_H
#define FAINT_SELECTION_QUERY_H

#include <string>
#include <utility>

namespace faint {

class SelectionQuery {
 public:
  SelectionQuery() = default;
  explicit SelectionQuery(std::string expression)
      : expression_(std::move(expression)) {}

  const std::string& str() const noexcept { return expression_; }
  bool empty() const noexcept { return expression_.empty(); }

 private:
  std::string expression_;
};

}  // namespace faint

#endif  // FAINT_SELECTION_QUERY_H

