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

  void update(size_t step, const double * state,
              dust::rng_state_t<real_t>& rng_state,
              double * state_next) {
    double S = state[0];
    double I = state[1];
    double R = state[2];
    double cumulative_incidence = state[3];

    double N = S + I + R;

    double p_SI = 1 - std::exp(-(data_.beta) * I / (double) N);
    double p_IR = 1 - std::exp(-(data_.gamma));
    double n_IR = dust::distr::rbinom(rng_state, I, p_IR * data_.dt);
    double n_SI = dust::distr::rbinom(rng_state, S, p_SI * data_.dt);

    state_next[0] = S - n_SI;
    state_next[1] = I + n_SI - n_IR;
    state_next[2] = R + n_IR;
    state_next[3] = cumulative_incidence + n_SI;
  }

private:
  init_t data_;
};

#include <cpp11/list.hpp>

inline double get_or_set_default(cpp11::list data, const char * name,
                                 double default_value) {
  SEXP value = data[name];
  return value == R_NilValue ? default_value : cpp11::as_cpp<double>(value);
}

template <>
sir::init_t dust_data<sir>(cpp11::list data) {
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
cpp11::sexp dust_info<sir>(const sir::init_t& data) {
  using namespace cpp11::literals;
  // Information about state order
  cpp11::writable::strings vars({"S", "I", "R", "inc"});
  // Information about parameter values
  cpp11::list pars = cpp11::writable::list({"beta"_nm = data.beta,
                                            "gamma"_nm = data.gamma});
  return cpp11::writable::list({"vars"_nm = vars, "pars"_nm = pars});
}
