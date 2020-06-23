models <- list(
  sir = dust::dust(mcstate_file("example/sir/dust_sir.cpp"), type = "sir",
                   quiet = TRUE))


example_sir <- function() {
  model <- models$sir
  sir <- model$new(data = NULL, step = 0,
                   n_particles = 1)
  y0 <- sir$state()

  compare <- function(state, prev_state, observed, pars = NULL) {
    if (is.null(pars)) {
      exp_noise = 1e6
    } else {
      exp_noise = pars['exp_noise']
    }
    incidence_modelled <- prev_state[1, ] - state[1, ]
    incidence_observed <- observed$incidence
    lambda <- incidence_modelled +
      rexp(n = length(incidence_modelled), rate = exp_noise)
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
  }

  inv_dt <- 4
  day <- seq(1, 100)
  incidence <- rep(NA, length(day))
  history <- array(NA_real_, c(3, 1, day[length(day)] + 1))
  history[, , 1] <- sir$state()

  for (i in day) {
    state_start <- sir$state()
    sir$run(i * inv_dt)
    state_end <- sir$state()
    history[, , i] <- state_end
    # Reduction in S
    incidence[i] <- state_start[1, 1] - state_end[1, 1]
  }

  data_raw <- data.frame(day = day, incidence = incidence)
  data <- particle_filter_data(data_raw, "day", 4)

  list(model = model, compare = compare, y0 = y0,
       data_raw = data_raw, data = data, history = history)
}


data_frame <- function(...) {
  data.frame(..., stringsAsFactors = FALSE)
}
