#include "rarexsec/fit/Fitter.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TH1.h"
#include "TH1D.h"

namespace rarexsec::internal::fit {

Fitter::Fitter(const std::string &signal_process_label) : signal_label_(signal_process_label) {}

Fitter::Fitter(const Fitter &o)
    : channels_(o.channels_),
      all_channels_(o.all_channels_),
      all_processes_(o.all_processes_),
      norm_nuis_(o.norm_nuis_),
      shape_nuis_(o.shape_nuis_),
      signal_label_(o.signal_label_),
      sigma_ref_pb_(o.sigma_ref_pb_),
      mu_lo_(o.mu_lo_),
      mu_hi_(o.mu_hi_),
      eps_(o.eps_),
      n_pars_(o.n_pars_),
      par_names_(o.par_names_),
      par_is_norm_(o.par_is_norm_) {}

Fitter &Fitter::operator=(const Fitter &o) {
  if (this == &o) return *this;
  channels_ = o.channels_;
  all_channels_ = o.all_channels_;
  all_processes_ = o.all_processes_;
  norm_nuis_ = o.norm_nuis_;
  shape_nuis_ = o.shape_nuis_;
  signal_label_ = o.signal_label_;
  sigma_ref_pb_ = o.sigma_ref_pb_;
  mu_lo_ = o.mu_lo_;
  mu_hi_ = o.mu_hi_;
  eps_ = o.eps_;
  n_pars_ = o.n_pars_;
  par_names_ = o.par_names_;
  par_is_norm_ = o.par_is_norm_;
  return *this;
}

Fitter::~Fitter() { clear_(); }


Fitter::Process::Process(const Process &o) : name(o.name), is_signal(o.is_signal) {
  if (o.nominal) {
    auto *cl = static_cast<TH1D *>(o.nominal->Clone());
    if (!cl) throw std::runtime_error("Process copy: clone failed");
    cl->SetDirectory(nullptr);
    nominal.reset(cl);
  }
}

Fitter::Process &Fitter::Process::operator=(const Process &o) {
  if (this == &o) return *this;
  name = o.name;
  is_signal = o.is_signal;
  if (o.nominal) {
    auto *cl = static_cast<TH1D *>(o.nominal->Clone());
    if (!cl) throw std::runtime_error("Process copy: clone failed");
    cl->SetDirectory(nullptr);
    nominal.reset(cl);
  } else {
    nominal.reset();
  }
  return *this;
}

Fitter::Channel::Channel(const Channel &o)
    : name(o.name), processes(o.processes), nbins(o.nbins) {
  if (o.data) {
    auto *cl = static_cast<TH1D *>(o.data->Clone());
    if (!cl) throw std::runtime_error("Channel copy: clone failed");
    cl->SetDirectory(nullptr);
    data.reset(cl);
  }
}

Fitter::Channel &Fitter::Channel::operator=(const Channel &o) {
  if (this == &o) return *this;
  name = o.name;
  nbins = o.nbins;
  processes = o.processes;
  if (o.data) {
    auto *cl = static_cast<TH1D *>(o.data->Clone());
    if (!cl) throw std::runtime_error("Channel copy: clone failed");
    cl->SetDirectory(nullptr);
    data.reset(cl);
  } else {
    data.reset();
  }
  return *this;
}

Fitter::ShapeNuisance::ShapeNuisance(const ShapeNuisance &o) : name(o.name), index(o.index) {
  for (const auto &kv : o.updown) {
    std::unique_ptr<TH1D> up;
    std::unique_ptr<TH1D> dn;
    if (kv.second.first) {
      auto *cl = static_cast<TH1D *>(kv.second.first->Clone());
      if (!cl) throw std::runtime_error("ShapeNuisance copy: clone failed (up)");
      cl->SetDirectory(nullptr);
      up.reset(cl);
    }
    if (kv.second.second) {
      auto *cl = static_cast<TH1D *>(kv.second.second->Clone());
      if (!cl) throw std::runtime_error("ShapeNuisance copy: clone failed (down)");
      cl->SetDirectory(nullptr);
      dn.reset(cl);
    }
    updown.emplace(kv.first, std::make_pair(std::move(up), std::move(dn)));
  }
}

Fitter::ShapeNuisance &Fitter::ShapeNuisance::operator=(const ShapeNuisance &o) {
  if (this == &o) return *this;
  name = o.name;
  index = o.index;
  updown.clear();
  for (const auto &kv : o.updown) {
    std::unique_ptr<TH1D> up;
    std::unique_ptr<TH1D> dn;
    if (kv.second.first) {
      auto *cl = static_cast<TH1D *>(kv.second.first->Clone());
      if (!cl) throw std::runtime_error("ShapeNuisance copy: clone failed (up)");
      cl->SetDirectory(nullptr);
      up.reset(cl);
    }
    if (kv.second.second) {
      auto *cl = static_cast<TH1D *>(kv.second.second->Clone());
      if (!cl) throw std::runtime_error("ShapeNuisance copy: clone failed (down)");
      cl->SetDirectory(nullptr);
      dn.reset(cl);
    }
    updown.emplace(kv.first, std::make_pair(std::move(up), std::move(dn)));
  }
  return *this;
}

void Fitter::set_sigma_ref(double sigma_ref_pb) { sigma_ref_pb_ = sigma_ref_pb; }

double Fitter::sigma_ref() const { return sigma_ref_pb_; }

void Fitter::set_mu_bounds(double lo, double hi) {
  mu_lo_ = lo;
  mu_hi_ = hi;
}

void Fitter::set_yield_floor(double eps) { eps_ = (eps > 0.0 ? eps : 1e-12); }

void Fitter::add_channel(const std::string &channel, const TH1 *h_data) {
  if (!h_data) throw std::invalid_argument("add_channel: data histogram is null");
  if (channels_.count(channel)) throw std::runtime_error("channel already exists: " + channel);
  Channel ch;
  ch.name = channel;
  ch.data.reset(clone_as_th1d_(h_data, channel + "__data"));
  ch.nbins = ch.data->GetNbinsX();
  channels_.emplace(channel, std::move(ch));
}

void Fitter::add_process(const std::string &channel, const std::string &process, const TH1 *h_nominal,
                         bool is_signal) {
  auto it = channels_.find(channel);
  if (it == channels_.end()) throw std::runtime_error("add_process: unknown channel " + channel);
  if (!h_nominal) throw std::invalid_argument("add_process: nominal histogram is null");
  if (it->second.processes.count(process))
    throw std::runtime_error("add_process: process already exists in channel: " + process);
  auto h = std::unique_ptr<TH1D>(clone_as_th1d_(h_nominal, channel + "__" + process + "__nom"));
  ensure_same_binning_(*it->second.data, *h, "add_process(" + channel + "," + process + ")");
  Process p;
  p.name = process;
  p.is_signal = (is_signal || process == signal_label_);
  p.nominal = std::move(h);
  it->second.processes.emplace(process, std::move(p));
  all_channels_.insert(channel);
  all_processes_.insert(process);
}

void Fitter::mark_signal_process(const std::string &process) {
  signal_label_ = process;
  for (auto &kv : channels_) {
    for (auto &pkv : kv.second.processes) {
      pkv.second.is_signal = (pkv.first == signal_label_);
    }
  }
}

void Fitter::add_norm_systematic(const std::string &name, bool log_normal) {
  if (norm_nuis_.count(name)) throw std::runtime_error("norm nuisance already exists: " + name);
  NormNuisance nn;
  nn.name = name;
  nn.log_normal = log_normal;
  norm_nuis_.emplace(name, std::move(nn));
}

void Fitter::set_norm_effect(const std::string &name, const std::string &channel, const std::string &process,
                             double frac) {
  auto it = norm_nuis_.find(name);
  if (it == norm_nuis_.end()) throw std::runtime_error("unknown norm nuisance: " + name);
  if (!has_proc_(channel, process))
    throw std::runtime_error("set_norm_effect: unknown (channel, process): " + channel + "," + process);
  if (frac < 0) throw std::invalid_argument("set_norm_effect: frac must be >= 0");
  it->second.frac[CPKey{channel, process}] = frac;
}

void Fitter::add_shape_systematic(const std::string &name) {
  if (shape_nuis_.count(name)) throw std::runtime_error("shape nuisance already exists: " + name);
  ShapeNuisance sn;
  sn.name = name;
  shape_nuis_.emplace(name, std::move(sn));
}

void Fitter::set_shape_effect(const std::string &name, const std::string &channel, const std::string &process,
                              const TH1 *h_up, const TH1 *h_down) {
  auto it = shape_nuis_.find(name);
  if (it == shape_nuis_.end()) throw std::runtime_error("unknown shape nuisance: " + name);
  if (!has_proc_(channel, process))
    throw std::runtime_error("set_shape_effect: unknown (channel, process): " + channel + "," + process);
  if (!h_up || !h_down) throw std::invalid_argument("set_shape_effect: h_up/h_down must be non-null");
  auto &nom = *channels_.at(channel).processes.at(process).nominal;
  auto up = std::unique_ptr<TH1D>(clone_as_th1d_(h_up, channel + "__" + process + "__" + name + "__up"));
  auto dn = std::unique_ptr<TH1D>(clone_as_th1d_(h_down, channel + "__" + process + "__" + name + "__down"));
  ensure_same_binning_(nom, *up, "set_shape_effect(up:" + channel + "," + process + "," + name + ")");
  ensure_same_binning_(nom, *dn, "set_shape_effect(down:" + channel + "," + process + "," + name + ")");
  it->second.updown[CPKey{channel, process}] = std::make_pair(std::move(up), std::move(dn));
}

Fitter::FitResult Fitter::fit(const std::string &minimizer, const std::string &algo, bool verbose) {
  if (channels_.empty()) throw std::runtime_error("fit: no channels added");
  if (!has_any_signal_()) throw std::runtime_error("fit: no signal process marked");
  build_parameter_indexing_();
  std::unique_ptr<ROOT::Math::Minimizer> min{ROOT::Math::Factory::CreateMinimizer(minimizer.c_str(),
                                                                                 algo.c_str())};
  if (!min) throw std::runtime_error("failed to create ROOT::Math::Minimizer");
  min->SetPrintLevel(verbose ? 1 : 0);
  min->SetStrategy(1);
  min->SetMaxFunctionCalls(100000);
  min->SetMaxIterations(100000);
  min->SetTolerance(1e-4);
  ROOT::Math::Functor f(this, &Fitter::nll_, n_pars_);
  min->SetFunction(f);
  min->SetLimitedVariable(0, "mu", std::clamp(guess_mu_(), mu_lo_, mu_hi_), 0.1, mu_lo_, mu_hi_);
  for (std::size_t i = 1; i < n_pars_; ++i)
    min->SetVariable(static_cast<int>(i), par_names_[i].c_str(), 0.0, 0.1);
  bool ok = min->Minimize();
  FitResult fr;
  fr.status = min->Status();
  fr.nll = min->MinValue() / 2.0;
  const double *xs = min->X();
  const double *xe = min->Errors();
  if (ok && xs) {
    fr.mu = xs[0];
    fr.mu_err_sym = (xe ? xe[0] : std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 1; i < n_pars_; ++i) {
      fr.nuis_values[par_names_[i]] = xs[i];
      if (xe) fr.nuis_errors[par_names_[i]] = xe[i];
    }
  }
  min->Hesse();
  if (min->Errors()) fr.mu_err_sym = min->Errors()[0];
  return fr;
}

std::vector<std::pair<double, double>> Fitter::scan_delta_nll(double mu_min, double mu_max, int npts,
                                                              const std::string &minimizer,
                                                              const std::string &algo,
                                                              bool verbose) {
  if (mu_min >= mu_max) throw std::invalid_argument("scan_delta_nll: mu_min < mu_max required");
  if (npts < 3) throw std::invalid_argument("scan_delta_nll: npts >= 3 required");
  build_parameter_indexing_();
  std::unique_ptr<ROOT::Math::Minimizer> min{ROOT::Math::Factory::CreateMinimizer(minimizer.c_str(),
                                                                                 algo.c_str())};
  if (!min) throw std::runtime_error("failed to create ROOT::Math::Minimizer");
  min->SetPrintLevel(verbose ? 1 : 0);
  min->SetStrategy(1);
  min->SetMaxFunctionCalls(200000);
  min->SetMaxIterations(200000);
  min->SetTolerance(1e-4);
  ROOT::Math::Functor f(this, &Fitter::nll_, n_pars_);
  min->SetFunction(f);
  double mu0 = std::clamp(guess_mu_(), mu_lo_, mu_hi_);
  min->SetLimitedVariable(0, "mu", mu0, 0.1, mu_lo_, mu_hi_);
  min->FixVariable(0);
  for (std::size_t i = 1; i < n_pars_; ++i)
    min->SetVariable(static_cast<int>(i), par_names_[i].c_str(), 0.0, 0.1);
  double nll_min_global = get_nll_min_free_mu_(minimizer, algo, verbose);
  std::vector<std::pair<double, double>> out;
  out.reserve(npts);
  for (int ip = 0; ip < npts; ++ip) {
    const double mu = mu_min + (mu_max - mu_min) * (double(ip) / double(npts - 1));
    min->SetVariableValue(0, mu);
    min->FixVariable(0);
    min->Minimize();
    double nll = min->MinValue() / 2.0;
    out.emplace_back(mu, std::max(0.0, nll - nll_min_global));
  }
  return out;
}

double Fitter::cross_section_pb(const FitResult &fr) const { return fr.mu * sigma_ref_pb_; }

double Fitter::cross_section_err_sym_pb(const FitResult &fr) const { return fr.mu_err_sym * sigma_ref_pb_; }

TH1D *Fitter::clone_as_th1d_(const TH1 *h, const std::string &new_name) {
  if (!h) return nullptr;
  TH1 *hc = dynamic_cast<TH1 *>(h->Clone(new_name.c_str()));
  if (!hc) throw std::runtime_error("clone failed for histogram " + std::string(h->GetName()));
  TH1D *hd = dynamic_cast<TH1D *>(hc);
  if (!hd) {
    TH1D *hnew = new TH1D(new_name.c_str(), h->GetTitle(), h->GetNbinsX(), h->GetXaxis()->GetXmin(),
                          h->GetXaxis()->GetXmax());
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      hnew->SetBinContent(i, h->GetBinContent(i));
      hnew->SetBinError(i, h->GetBinError(i));
    }
    delete hc;
    hd = hnew;
  }
  hd->SetDirectory(nullptr);
  return hd;
}

void Fitter::ensure_same_binning_(const TH1 &a, const TH1 &b, const std::string &ctx) {
  if (a.GetNbinsX() != b.GetNbinsX() || a.GetXaxis()->GetXmin() != b.GetXaxis()->GetXmin() ||
      a.GetXaxis()->GetXmax() != b.GetXaxis()->GetXmax()) {
    throw std::runtime_error("binning mismatch in " + ctx);
  }
}

bool Fitter::has_proc_(const std::string &ch, const std::string &pr) const {
  auto it = channels_.find(ch);
  if (it == channels_.end()) return false;
  return it->second.processes.count(pr) != 0;
}

bool Fitter::has_any_signal_() const {
  for (auto const &ckv : channels_)
    for (auto const &pkv : ckv.second.processes)
      if (pkv.second.is_signal) return true;
  return false;
}

void Fitter::clear_() {
  channels_.clear();
  norm_nuis_.clear();
  shape_nuis_.clear();
  all_channels_.clear();
  all_processes_.clear();
  par_names_.clear();
  par_is_norm_.clear();
  n_pars_ = 0;
}

void Fitter::build_parameter_indexing_() {
  par_names_.clear();
  par_is_norm_.clear();
  par_names_.push_back("mu");
  par_is_norm_.push_back(0);
  for (auto &kv : norm_nuis_) {
    kv.second.index = static_cast<int>(par_names_.size());
    par_names_.push_back("theta_norm_" + kv.first);
    par_is_norm_.push_back(1);
  }
  for (auto &kv : shape_nuis_) {
    kv.second.index = static_cast<int>(par_names_.size());
    par_names_.push_back("theta_shape_" + kv.first);
    par_is_norm_.push_back(2);
  }
  n_pars_ = par_names_.size();
}

double Fitter::guess_mu_() const {
  double s = 0.0, b = 0.0, d = 0.0;
  for (auto const &ckv : channels_) {
    const auto &ch = ckv.second;
    d += ch.data->Integral(1, ch.data->GetNbinsX());
    for (auto const &pkv : ch.processes) {
      const auto &p = pkv.second;
      const double y = p.nominal->Integral(1, p.nominal->GetNbinsX());
      if (p.is_signal)
        s += y;
      else
        b += y;
    }
  }
  double mu = (s > 0.0 ? (d - b) / s : 1.0);
  if (!std::isfinite(mu)) mu = 1.0;
  return std::clamp(mu, mu_lo_, mu_hi_);
}

double Fitter::nll_(const double *x) const {
  const double mu = std::clamp(x[0], mu_lo_, mu_hi_);
  double logl = 0.0;
  for (auto const &kv : norm_nuis_) {
    const double th = x[kv.second.index];
    logl += -0.5 * th * th;
  }
  for (auto const &kv : shape_nuis_) {
    const double th = x[kv.second.index];
    logl += -0.5 * th * th;
  }
  for (auto const &ckv : channels_) {
    const Channel &ch = ckv.second;
    for (int ib = 1; ib <= ch.nbins; ++ib) {
      double nu = 0.0;
      for (auto const &pkv : ch.processes) {
        const Process &proc = pkv.second;
        double y = proc.nominal->GetBinContent(ib);
        for (auto const &snkv : shape_nuis_) {
          const ShapeNuisance &sn = snkv.second;
          auto it = sn.updown.find(CPKey{ch.name, proc.name});
          if (it != sn.updown.end()) {
            const double th = x[sn.index];
            const double yup = it->second.first->GetBinContent(ib);
            const double ydn = it->second.second->GetBinContent(ib);
            const double delta = 0.5 * (yup - ydn);
            y += th * delta;
          }
        }
        if (y < 0.0) y = 0.0;
        double scale = 1.0;
        for (auto const &nnkv : norm_nuis_) {
          const NormNuisance &nn = nnkv.second;
          auto it = nn.frac.find(CPKey{ch.name, proc.name});
          if (it != nn.frac.end()) {
            const double th = x[nn.index];
            const double f = it->second;
            if (nn.log_normal)
              scale *= std::exp(std::log(1.0 + f) * th);
            else
              scale *= std::max(0.0, 1.0 + f * th);
          }
        }
        const double term = (proc.is_signal ? mu * scale * y : scale * y);
        nu += term;
      }
      const double nobs = ch.data->GetBinContent(ib);
      const double ex = (nu > eps_ ? nu : eps_);
      if (nobs > 0.0)
        logl += nobs * std::log(ex) - ex;
      else
        logl += -ex;
    }
  }
  return -2.0 * logl;
}

double Fitter::get_nll_min_free_mu_(const std::string &minimizer, const std::string &algo, bool verbose) {
  std::unique_ptr<ROOT::Math::Minimizer> min{ROOT::Math::Factory::CreateMinimizer(minimizer.c_str(),
                                                                                 algo.c_str())};
  min->SetPrintLevel(verbose ? 1 : 0);
  min->SetStrategy(1);
  min->SetMaxFunctionCalls(100000);
  min->SetMaxIterations(100000);
  min->SetTolerance(1e-4);
  ROOT::Math::Functor f(this, &Fitter::nll_, n_pars_);
  min->SetFunction(f);
  min->SetLimitedVariable(0, "mu", std::clamp(guess_mu_(), mu_lo_, mu_hi_), 0.1, mu_lo_, mu_hi_);
  for (std::size_t i = 1; i < n_pars_; ++i)
    min->SetVariable(static_cast<int>(i), par_names_[i].c_str(), 0.0, 0.1);
  min->Minimize();
  return min->MinValue() / 2.0;
}

} // namespace rarexsec::internal::fit

