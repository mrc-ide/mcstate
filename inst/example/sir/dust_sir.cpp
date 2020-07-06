// NB: includes such as dust/rng.hpp are automatic

class sir {
public:
  typedef int int_t;
  typedef double real_t;
  struct init_t {
    double S0;
    double I0;
    double R0;
    double beta;
    double gamma;
    double dt;
  };

  sir(const init_t& data) : data_(data) {
  }

  size_t size() const {
    return 4;
  }

  std::vector<double> initial(size_t step) {
    std::vector<double> ret = {data_.S0, data_.I0, data_.R0, 0};
    return ret;
  }

  void update(size_t step, const std::vector<double> state,
              dust::RNG<double, int>& rng,
              std::vector<double>& state_next) {
    double S = state[0];
    double I = state[1];
    double R = state[2];
    double cumulative_incidence = state[3];

    double N = S + I + R;

    double p_SI = 1 - std::exp(-(data_.beta) * I / (double) N);
    double p_IR = 1 - std::exp(-(data_.gamma));
    double n_IR = rng.rbinom(round(I), p_IR * data_.dt);
    double n_SI = rng.rbinom(round(S), p_SI * data_.dt);

    state_next[0] = S - n_SI;
    state_next[1] = I + n_SI - n_IR;
    state_next[2] = R + n_IR;
    state_next[3] = cumulative_incidence + n_SI;
  }

private:
  init_t data_;
};

#include <Rcpp.h>

inline double get_or_set_default(Rcpp::List& data, const char * name,
                                 double default_value) {
  double value = default_value;
  if (data.containsElementNamed(name)) {
    value = data[name];
  }
  return value;
}

template <>
sir::init_t dust_data<sir>(Rcpp::List data) {
  // Initial state values
  double S0 = get_or_set_default(data, "S0", 1000.0);
  double I0 = get_or_set_default(data, "I0", 10.0);
  double R0 = get_or_set_default(data, "R0", 0.0);
  double beta = get_or_set_default(data, "beta", 0.2);
  double gamma = get_or_set_default(data, "gamma", 0.1);
  double dt = 0.25;

  return sir::init_t{S0, I0, R0, beta, gamma, dt};
}

template <>
Rcpp::RObject dust_info<sir>(const sir::init_t& data) {
  // Information about state order
  Rcpp::CharacterVector vars = Rcpp::CharacterVector::create("S", "I", "R");
  // Information about parameter values
  Rcpp::List pars = Rcpp::List::create(Rcpp::Named("beta") = data.beta,
                                       Rcpp::Named("gamma") = data.gamma);
  return Rcpp::List::create(Rcpp::Named("vars") = vars,
                            Rcpp::Named("pars") = pars);
}
