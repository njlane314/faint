#pragma once

#include <TH1D.h>
#include <TMatrixDSym.h>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace rarexsec { struct Entry; }

namespace rarexsec::systpack {

struct Config {
  bool include_ext = false;
  bool use_ppfx = true;
  bool use_genie = true;
  bool use_reint = true;
  int N_ppfx = 600;
  int N_genie = 500;
  int N_reint = 100;
  std::string ppfx_branch = "weightsPPFX";
  std::string ppfx_cv_branch = "ppfx_cv";
  std::string genie_branch = "weightsGenie";
  std::string genie_cv_branch = "weightSplineTimesTune";
  std::string reint_branch = "weightsReint";
  double ushort_scale = 1.0 / 1000.0;
  std::string value_col = "x";
  std::string weight_col = "w_nominal";
};

struct Result {
  std::unique_ptr<TH1D> H_pred;
  std::map<std::string, TMatrixDSym> sources;
  TMatrixDSym total;
};

class SystematicsPack {
public:
  explicit SystematicsPack(Config cfg);
  Result build(const TH1D& model,
               const std::vector<const rarexsec::Entry*>& mc_entries,
               const std::vector<const rarexsec::Entry*>& ext_entries) const;
private:
  Config cfg_;
};

} // namespace rarexsec::systpack
