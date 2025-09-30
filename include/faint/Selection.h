#ifndef FAINT_SELECTION_H
#define FAINT_SELECTION_H

#include <string>

namespace faint {

class Selection {
 public:
  Selection();
  explicit Selection(std::string expression);

  const std::string& str() const noexcept;
  bool empty() const noexcept;

 private:
  std::string expression_;
};

}  // namespace faint

#endif  // FAINT_SELECTION_H
