models <- list(
  sir = dust::dust(mcstate_file("example/sir/dust_sir.cpp"), type = "sir",
                   quiet = TRUE))


example_sir <- function() {
  model <- models$sir
  sir <- model$new(data = NULL, step = 0, n_particles = 1)
  y0 <- sir$state()

  compare <- function(state, prev_state, observed, pars = NULL) {
    if (is.null(observed$incidence)) {
      return(NULL)
    }
    exp_noise <- pars[["exp_noise"]] %||% 1e6
    ## This is on the *filtered* state (i.e., returned by run())
    incidence_modelled <-
      state[1, , drop = TRUE] - prev_state[1, , drop = TRUE]
    incidence_observed <- observed$incidence
    lambda <- incidence_modelled +
      rexp(n = length(incidence_modelled), rate = exp_noise)
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
  }
  inv_dt <- 4
  day <- seq(1, 100)
  incidence <- rep(NA, length(day))
  history <- array(NA_real_, c(4, 1, length(day) + 1))
  history[, , 1] <- sir$state()

  for (i in day) {
    state_start <- sir$state()
    state_end <- sir$run(i * inv_dt)
    history[, , i + 1] <- state_end
    incidence[i] <- state_end[4, 1] - state_start[4, 1]
  }

  data_raw <- data.frame(day = day, incidence = incidence)
  data <- particle_filter_data(data_raw, "day", 4)
  index <- function(info) {
    list(run = 4L, state = 1:3)
  }

  list(model = model, compare = compare, y0 = y0,
       data_raw = data_raw, data = data, history = history,
       index = index)
}


data_frame <- function(...) {
  data.frame(..., stringsAsFactors = FALSE)
}
