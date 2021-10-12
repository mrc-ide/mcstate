context("deterministic")

test_that("Can run the deterministic filter", {
  dat <- example_sir()
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  set.seed(1)
  ll <- p$run(dat$pars$model(dat$pars$initial()))
  expect_equal(ll, -254.317276732767)
})


test_that("Can control starting point of simulation", {
  initial <- function(info, n_particles, pars) {
    list(step = pars$initial)
  }

  dat <- example_sir()

  data <- dat$data
  offset <- 400
  data[c("step_start", "step_end")] <-
    data[c("step_start", "step_end")] + offset
  data$step_start[[1]] <- 0

  pars <- dat$pars$model(dat$pars$initial())

  ## The usual version:
  p1 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                   index = dat$index)
  set.seed(1)
  ll1 <- p1$run(pars)

  ## Tuning the start date
  p2 <- particle_deterministic$new(data, dat$model, dat$compare,
                                   index = dat$index, initial = initial)
  set.seed(1)
  ll2 <- p2$run(c(pars, list(initial = as.integer(offset))))
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


test_that("can't restart deterministic filter", {
  dat <- example_sir()
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  expect_error(
    p$run(list(), save_restart = c(100, 200)),
    "'save_restart' cannot be used with particle_deterministic")
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
  expect_equal(p$history(1L), array_drop(h, 2))
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
  data_nested <- structure(dat$data,
                           class = c("particle_filter_data",
                                     "particle_filter_data_nested"))
  expect_error(
    particle_deterministic$new(data_nested, dat$model, dat$compare, dat$index),
    "nested mode not yet supported")
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


test_that("initial handles corner cases", {
  initial1 <- function(info, n_particles, pars) {
    rep(pars, 4)
  }
  initial2 <- function(info, n_particles, pars) {
    list(state = rep(pars, 4))
  }
  initial3 <- function(info, n_particles, pars) {
    list(state = rep(pars, 4), step = pars * 2)
  }

  pars <- list(1, 2, 3)
  info <- vector("list", length(pars))
  p1 <- deterministic_initial(pars, initial1, info)
  p2 <- deterministic_initial(pars, initial2, info)
  p3 <- deterministic_initial(pars, initial3, info)

  m <- matrix(rep(1:3, each = 4), 4, 3)
  expect_equal(p1, list(state = m))
  expect_equal(p2, p1)
  expect_equal(p3, list(state = m, step = c(2, 4, 6)))
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


test_that("initial passes args as expected, for multipar case", {
  dat <- example_sir()
  initial <- mockery::mock(c(1000, 10, 0, 0, 0), cycle = TRUE)
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                  index = dat$index, initial = initial)
  pars <- list(list(beta = 0.21, gamma = 0.11),
               list(beta = 0.22, gamma = 0.12),
               list(beta = 0.23, gamma = 0.13))
  ll <- p$run_many(pars)
  info <- lapply(pars, function(p) dat$model$new(p, 0, 1)$info())
  mockery::expect_called(initial, 3L)
  expect_equal(mockery::mock_args(initial)[[1L]],
               list(info[[1L]], 1L, pars[[1L]]))
  expect_equal(mockery::mock_args(initial)[[2L]],
               list(info[[2L]], 1L, pars[[2L]]))
  expect_equal(mockery::mock_args(initial)[[3L]],
               list(info[[3L]], 1L, pars[[3L]]))
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


test_that("initial passes args as expected, for multipar case", {
  ## As normal, but we require the exp_noise parameter
  compare <- function(state, observed, pars = NULL) {
    if (is.na(observed$incidence)) {
      return(NULL)
    }
    exp_noise <- pars$exp_noise
    incidence_modelled <- state[1L, , drop = TRUE]
    incidence_observed <- observed$incidence
    lambda <- incidence_modelled +
      rexp(n = length(incidence_modelled), rate = exp_noise)
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
  }

  dat <- example_sir()
  p <- particle_deterministic$new(dat$data, dat$model, compare,
                                  index = dat$index, initial = dat$initial)
  pars <- list(list(exp_noise = Inf),
               list(exp_noise = 1e6))
  ll <- p$run_many(pars)
  expect_equal(ll[[1]], ll[[2]], tolerance = 1e-5)
  expect_false(identical(ll[[1]], ll[[2]]))
})


test_that("reconstruct deterministic filter from inputs", {
  dat <- example_sir()
  p1 <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  inputs <- p1$inputs()
  p2 <- particle_filter_from_inputs(inputs)
  expect_s3_class(p2, class(p1))
  expect_equal(p1$run(), p2$run())
})


test_that("Can run parallel mcmc with deterministic model", {
  dat <- example_sir()
  n_steps <- 30L
  n_chains <- 3L
  control <- pmcmc_control(n_steps, save_trajectories = TRUE,
                           n_workers = 2L, n_chains = n_chains,
                           save_state = TRUE)
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare, dat$index)
  res <- pmcmc(dat$pars, p, NULL, control)
  expect_s3_class(res, "mcstate_pmcmc")
  expect_equal(nrow(res$pars), n_chains * (n_steps + 1))
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
