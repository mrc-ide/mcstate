test_that("Trivial multistage deterministic same as simple", {
  dat <- example_sir()
  index <- function(info) {
    list(run = 5L, state = c(S = 1, I = 2, R = 3))
  }

  pars_base <- list(beta = 0.2, gamma = 0.1, compare = list(exp_noise = Inf))
  pars <- multistage_parameters(pars_base, list())

  set.seed(1)
  filter1 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                        index = index)
  ll1 <- filter1$run(pars_base, save_history = TRUE)

  set.seed(2)
  filter2 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                        index = index)
  ll2 <- filter2$run(pars, save_history = TRUE)

  expect_identical(ll1, ll2)

  ## expect_identical(
  ##   r6_private(filter1)$last_model$rng_state(),
  ##   r6_private(filter2)$last_model$rng_state())

  ## Not totally sure how we don't have rownames here.
  expect_identical(filter1$history(), filter2$history())
})


test_that("An effectless multistage filter is identical to single stage", {
  dat <- example_sir()
  index <- function(info) {
    list(run = 5L, state = c(S = 1, I = 2, R = 3))
  }

  epochs <- list(multistage_epoch(10), multistage_epoch(20))
  pars_base <- list(beta = 0.2, gamma = 0.1, compare = list(exp_noise = Inf))
  pars <- multistage_parameters(pars_base, epochs)
  restart <- c(5, 15, 25)

  set.seed(1)
  filter1 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                        index = index)
  ll1 <- filter1$run(pars_base, save_history = TRUE, save_restart = restart)

  set.seed(1)
  filter2 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                        index = index)
  ll2 <- filter2$run(pars, save_history = TRUE, save_restart = restart)

  ## To get these to be identical we must loop over the data in the
  ## other order (by pars-within-time, then time); see
  ## test-deterministic.R for a similar problem (but note that the
  ## actual states here are spot on)
  expect_equal(ll1, ll2, tolerance = 1e-11)
  ## expect_identical(
  ##   r6_private(filter1)$last_model$rng_state(),
  ##   r6_private(filter2)$last_model$rng_state())

  expect_identical(filter1$history(), filter2$history())
  expect_identical(filter1$restart_state(), filter2$restart_state())
})
