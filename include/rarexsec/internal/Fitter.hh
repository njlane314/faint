#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

class TH1;
class TH1D;

namespace rarexsec::internal::fit {

class Fitter {
public:
  struct FitResult {
    int status = -1;
    double nll = std::numeric_limits<double>::quiet_NaN();
    double mu = std::numeric_limits<double>::quiet_NaN();
    double mu_err_sym = std::numeric_limits<double>::quiet_NaN();
    double mu_err_lo = std::numeric_limits<double>::quiet_NaN();
    double mu_err_hi = std::numeric_limits<double>::quiet_NaN();
    std::map<std::string, double> nuis_values;
    std::map<std::string, double> nuis_errors;
  };

  struct CPKey {
    std::string ch;
    std::string pr;
    bool operator<(const CPKey &o) const { return std::tie(ch, pr) < std::tie(o.ch, o.pr); }
  };

  explicit Fitter(const std::string &signal_process_label = "signal");
  ~Fitter();

  void set_sigma_ref(double sigma_ref_pb);
  double sigma_ref() const;
  void set_mu_bounds(double lo, double hi);
  void set_yield_floor(double eps);

  void add_channel(const std::string &channel, const TH1 *h_data);
  void add_process(const std::string &channel, const std::string &process, const TH1 *h_nominal,
                   bool is_signal = false);
  void mark_signal_process(const std::string &process);

  void add_norm_systematic(const std::string &name, bool log_normal = true);
  void set_norm_effect(const std::string &name, const std::string &channel, const std::string &process,
                       double frac);

  void add_shape_systematic(const std::string &name);
  void set_shape_effect(const std::string &name, const std::string &channel, const std::string &process,
                        const TH1 *h_up, const TH1 *h_down);

  FitResult fit(const std::string &minimizer = "Minuit2", const std::string &algo = "Migrad", bool verbose = false);

  std::vector<std::pair<double, double>> scan_delta_nll(double mu_min, double mu_max, int npts,
                                                        const std::string &minimizer = "Minuit2",
                                                        const std::string &algo = "Migrad",
                                                        bool verbose = false);

  double cross_section_pb(const FitResult &fr) const;
  double cross_section_err_sym_pb(const FitResult &fr) const;

private:
  struct Process {
    std::string name;
    bool is_signal = false;
    std::unique_ptr<TH1D> nominal;
  };

  struct Channel {
    std::string name;
    std::unique_ptr<TH1D> data;
    std::map<std::string, Process> processes;
    int nbins = 0;
  };

  struct NormNuisance {
    std::string name;
    bool log_normal = true;
    std::map<CPKey, double> frac;
    int index = -1;
  };

  struct ShapeNuisance {
    std::string name;
    std::map<CPKey, std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>>> updown;
    int index = -1;
  };

  static TH1D *clone_as_th1d_(const TH1 *h, const std::string &new_name);
  static void ensure_same_binning_(const TH1 &a, const TH1 &b, const std::string &ctx);

  bool has_proc_(const std::string &ch, const std::string &pr) const;
  bool has_any_signal_() const;
  void clear_();
  void build_parameter_indexing_();
  double guess_mu_() const;
  double nll_(const double *x) const;
  double get_nll_min_free_mu_(const std::string &minimizer, const std::string &algo, bool verbose);

  std::map<std::string, Channel> channels_;
  std::set<std::string> all_channels_;
  std::set<std::string> all_processes_;
  std::map<std::string, NormNuisance> norm_nuis_;
  std::map<std::string, ShapeNuisance> shape_nuis_;
  std::string signal_label_;
  double sigma_ref_pb_ = 1.0;
  double mu_lo_ = 0.0;
  double mu_hi_ = 10.0;
  double eps_ = 1e-9;
  std::size_t n_pars_ = 0;
  std::vector<std::string> par_names_;
  std::vector<int> par_is_norm_;
};

} // namespace rarexsec::internal::fit

