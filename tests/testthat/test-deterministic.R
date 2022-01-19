context("deterministic")

test_that("Can run the deterministic filter", {
  dat <- example_sir()
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  set.seed(1)
  ll <- p$run(dat$pars$model(dat$pars$initial()))
  expect_equal(ll, -254.317276732767)
})


test_that("Can control starting point of simulation", {
  dat <- example_sir()

  pars <- dat$pars$model(dat$pars$initial())

  ## The usual version:
  p1 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                   index = dat$index)
  set.seed(1)
  ll1 <- p1$run(pars)

  ## Tuning the start date
  data_raw <- dat$data_raw
  data_raw$day <- data_raw$day + 100
  data <- particle_filter_data(data_raw, "day", 4, 100)

  ## Tuning the start date
  p2 <- particle_deterministic$new(data, dat$model, dat$compare,
                                   index = dat$index)
  set.seed(1)
  ll2 <- p2$run(pars)
  expect_identical(ll1, ll2)
})


test_that("Control the initial conditions", {
  dat <- example_sir()

  initial <- function(info, n_particles, pars) {
    c(1000, pars$I0, 0, 0, 0)
  }
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                  index = dat$index, initial = initial)

  ll1 <- p$run(list(I0 = 200), save_history = TRUE)
  expect_equal(p$history()[, , 1],
               c(1000, 200, 0))
  ll2 <- p$run(list(I0 = 1), save_history = TRUE)
  expect_equal(p$history()[, , 1],
               c(1000, 1, 0))
  ll3 <- p$run(list(I0 = 10), save_history = TRUE)
  expect_equal(p$history()[, , 1],
               c(1000, 10, 0))

  expect_true(ll1 < ll3)
  expect_true(ll2 < ll3)
})


test_that("can collect restart from deterministic filter", {
  dat <- example_sir()
  set.seed(1)
  p1 <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  expect_error(p1$restart_state(),
               "Model has not yet been run")
  res1 <- p1$run(list())

  set.seed(1)
  p2 <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  res2 <- p2$run(list(), save_restart = c(25, 50), save_history = TRUE)

  expect_equal(res1, res2)
  expect_error(p1$restart_state(),
               "Can't get history as model was run with save_restart = NULL")

  h <- p2$history()
  r <- p2$restart_state()
  expect_equal(dim(r), c(5, 1, 2))
  expect_equal(h[, , 26], r[1:3, 1, 1])
  expect_equal(h[, , 51], r[1:3, 1, 2])

  ## Can resample correctly:
  expect_equal(p2$restart_state(1L), p2$restart_state())
})


test_that("Particle index control in filter is very limited", {
  dat <- example_sir()
  set.seed(1)
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  p$run(list(), save_history = TRUE, save_restart = 50)

  expect_equal(dim(p$history()), c(3, 1, 101))
  expect_equal(dim(p$history(1)), c(3, 1, 101))
  expect_error(
    dim(p$history(integer(0))),
    "Invalid value for 'index_particle' may only be 1 (or NULL)",
    fixed = TRUE)
  expect_error(
    dim(p$history(2)),
    "Invalid value for 'index_particle' may only be 1 (or NULL)",
    fixed = TRUE)
  expect_error(
    dim(p$history(c(1, 1))),
    "Invalid value for 'index_particle' may only be 1 (or NULL)",
    fixed = TRUE)

  ## This follows the behaviour of the particle filter, which is
  ## pretty weird really.
  expect_equal(dim(p$restart_state()), c(5, 1, 1))
  expect_equal(dim(p$restart_state(1)), c(5, 1, 1))
  expect_error(
    dim(p$restart_state(integer(0))),
    "Invalid value for 'index_particle' may only be 1 (or NULL)",
    fixed = TRUE)
  expect_error(
    dim(p$restart_state(2)),
    "Invalid value for 'index_particle' may only be 1 (or NULL)",
    fixed = TRUE)
  expect_error(
    dim(p$restart_state(c(1, 1))),
    "Invalid value for 'index_particle' may only be 1 (or NULL)",
    fixed = TRUE)
})


test_that("can't use min_log_likelihood with deterministic filter", {
  dat <- example_sir()
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  expect_error(
    p$run(list(), min_log_likelihood = -200),
    "'min_log_likelihood' cannot be used with particle_deterministic")
})


test_that("can't return state until filter run", {
  dat <- example_sir()
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  expect_error(p$state(),
               "Model has not yet been run")
})


test_that("can't return history until filter run, with save_history", {
  dat <- example_sir()
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  expect_error(p$history(),
               "Model has not yet been run")
  p$run()
  expect_error(p$history(),
               "Can't get history as model was run with save_history = FALSE")
})


test_that("extract history from deterministic filter", {
  dat <- example_sir()
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  p$run(save_history = TRUE)
  h <- p$history()
  expect_equal(dim(h), c(3, 1, 101)) # state, particles, time
  expect_equal(p$history(1L), h)
})


test_that("Validate inputs to the deterministic filter", {
  dat <- example_sir()
  expect_error(
    particle_deterministic$new(dat$data, NULL, dat$compare, dat$index),
    "'model' must be a dust_generator")
  expect_error(
    particle_deterministic$new(dat$data, dat$model, dat$compare, TRUE),
    "'index' must be function if not NULL")
  expect_error(
    particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index,
                               initial = 1),
    "'initial' must be function if not NULL")
})


test_that("Can run mcmc with deterministic filter", {
  dat <- example_sir()
  control <- pmcmc_control(100, save_trajectories = TRUE)
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  res <- pmcmc(dat$pars, p, NULL, control)
  expect_s3_class(res, "mcstate_pmcmc")
})


test_that("Can run deterministic filter without index", {
  dat <- example_sir()
  n_particles <- 42
  p1 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                   index = dat$index)

  compare2 <- function(state, ...) {
    dat$compare(state[5, , drop = FALSE], ...)
  }

  p2 <- particle_deterministic$new(dat$data, dat$model, compare2)

  set.seed(1)
  ll1 <- p1$run(save_history = TRUE)
  set.seed(1)
  ll2 <- p2$run(save_history = TRUE)
  expect_identical(ll1, ll2)

  expect_equal(dim(p1$history()), c(3, 1, 101))
  expect_equal(dim(p2$history()), c(5, 1, 101))
})


test_that("initial passes args as expected", {
  dat <- example_sir()
  initial <- mockery::mock(c(1000, 10, 0, 0, 0), cycle = TRUE)
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                  index = dat$index, initial = initial)
  pars <- list(beta = 0.21, gamma = 0.11)
  ll <- p$run(pars)
  info <- dat$model$new(pars, 0, 1)$info()
  mockery::expect_called(initial, 1L)
  expect_equal(mockery::mock_args(initial)[[1L]],
               list(info, 1L, pars))
})


test_that("Can reference by name in compare", {
  dat <- example_sir()
  p1 <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  set.seed(1)
  ll1 <- p1$run(dat$pars$model(dat$pars$initial()),
                save_history = TRUE)

  compare <- function(state, observed, pars = NULL) {
    if (is.na(observed$incidence)) {
      return(NULL)
    }
    if (is.null(pars$compare$exp_noise)) {
      exp_noise <- 1e6
    } else {
      exp_noise <- pars$compare$exp_noise
    }
    ## This is on the *filtered* state (i.e., returned by run())
    incidence_modelled <- state["incidence", , drop = TRUE]
    incidence_observed <- observed$incidence
    lambda <- incidence_modelled +
      rexp(n = length(incidence_modelled), rate = exp_noise)
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
  }

  index <- function(info) {
    list(run = c("S" = 1L, "incidence" = 5L),
         state = c("S" = 1L, "I" = 2L, "R" = 3L))
  }

  p2 <- particle_deterministic$new(dat$data, dat$model, compare, index)
  set.seed(1)
  ll2 <- p2$run(dat$pars$model(dat$pars$initial()),
                save_history = TRUE)
  expect_identical(ll2, ll1)
  h1 <- p1$history()
  h2 <- p2$history()

  expect_equal(h1, unname(h2))
  expect_equal(dimnames(h2), list(c("S", "I", "R"), NULL, NULL))
})


test_that("reconstruct deterministic filter from inputs", {
  dat <- example_sir()
  p1 <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  inputs <- p1$inputs()
  p2 <- particle_filter_from_inputs(inputs)
  expect_s3_class(p2, class(p1))
  expect_equal(p1$run(), p2$run(), tolerance = 1e-5)
})


test_that("Can run parallel mcmc with deterministic model", {
  dat <- example_sir()
  n_steps <- 30L
  n_chains <- 3L
  control <- pmcmc_control(n_steps, save_trajectories = TRUE,
                           n_workers = 2L, n_chains = n_chains,
                           n_threads_total = 2, save_state = TRUE)
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  res <- pmcmc(dat$pars, p, NULL, control)
  expect_s3_class(res, "mcstate_pmcmc")
  expect_equal(nrow(res$pars), n_chains * n_steps)
  expect_s3_class(res$predict$filter$model, "dust_generator")
})


test_that("Can change the number of threads", {
  dat <- example_sir()
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  expect_equal(p$set_n_threads(2L), 1L)
  expect_equal(p$set_n_threads(1L), 2L)
  p$run()
  expect_equal(p$set_n_threads(2L), 1L)
  expect_equal(p$set_n_threads(1L), 2L)
})


test_that("Cannot use previous initial condition approach", {
  initial <- function(info, n_particles, pars) {
    list(step = 2)
  }
  dat <- example_sir()
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                  index = dat$index, initial = initial)
  expect_error(p$run(), "Setting 'step' from initial no longer supported")
})


test_that("Can partially run a deterministic particle and resume", {
  dat <- example_sir()
  pars <- dat$pars$model(dat$pars$initial())

  set.seed(1)
  p1 <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  ll1 <- p1$run(pars, save_history = TRUE)

  set.seed(1)
  p2 <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  s2 <- p2$run_begin(list(pars), save_history = TRUE)
  ll2 <- s2$run()
  expect_identical(ll2, ll1)
  expect_identical(s2$history$value, p1$history())

  set.seed(1)
  p3 <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  s3 <- p3$run_begin(list(pars), save_history = TRUE)
  ll3_1 <- s3$step(50)
  ll3_2 <- s3$step(100)
  ## To get these to be identical we must loop over the data in the
  ## other order (by pars-within-time, then time)
  expect_equal(ll3_2, ll1, tolerance = 1e-11)
  expect_equal(s3$history$value, p1$history(), tolerance = 1e-11)
})


test_that("Can offset the initial likelihood", {
  dat <- example_sir()
  n_particles <- 42

  constant_ll <- function(pars) {
    10
  }

  set.seed(1)
  p1 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                   index = dat$index)
  ll1 <- p1$run(save_history = TRUE)

  set.seed(1)
  p2 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                            constant_log_likelihood = constant_ll,
                            index = dat$index)
  ll2 <- p2$run(save_history = TRUE)
  expect_equal(ll2, ll1 + 10)
  expect_identical(p1$history(), p2$history())
  expect_identical(p1$state(), p2$state())
})
