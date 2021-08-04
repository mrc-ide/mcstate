context("deterministic")

test_that("Can run the deterministic filter", {
  dat <- example_sir()
  p <- particle_nofilter$new(dat$data, dat$model, dat$compare, dat$index)
  set.seed(1)
  ll <- p$run(dat$pars$model(dat$pars$initial()))
  expect_equal(ll, -245.127512965178)
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
  p1 <- particle_nofilter$new(dat$data, dat$model, dat$compare,
                              index = dat$index)
  set.seed(1)
  ll1 <- p1$run(pars)

  ## Tuning the start date
  p2 <- particle_nofilter$new(data, dat$model, dat$compare,
                              index = dat$index, initial = initial)
  set.seed(1)
  ll2 <- p2$run(c(pars, list(initial = as.integer(offset))))
  expect_identical(ll1, ll2)
})


test_that("can't restart deterministic filter", {
  dat <- example_sir()
  p <- particle_nofilter$new(dat$data, dat$model, dat$compare, dat$index)
  expect_error(
    p$run(pars, save_restart = c(100, 200)),
    "'save_restart' cannot be used with the deterministic nofilter")
})


test_that("can't return state until filter run", {
  dat <- example_sir()
  p <- particle_nofilter$new(dat$data, dat$model, dat$compare, dat$index)
  expect_error(p$state(),
               "Model has not yet been run")
})


test_that("can't return history until filter run, with save_history", {
  dat <- example_sir()
  p <- particle_nofilter$new(dat$data, dat$model, dat$compare, dat$index)
  expect_error(p$history(),
               "Model has not yet been run")
  p$run()
  expect_error(p$history(),
               "Can't get history as model was run with save_history = FALSE")
})


test_that("extract history from deterministic filter", {
  dat <- example_sir()
  p <- particle_nofilter$new(dat$data, dat$model, dat$compare, dat$index)
  p$run(save_history = TRUE)
  h <- p$history()
  expect_equal(dim(h), c(3, 1, 101)) # state, particles, time
  expect_equal(p$history(1L), array_drop(h, 2))
})


test_that("Validate inputs to the deterministic filter", {
  dat <- example_sir()
  expect_error(
    particle_nofilter$new(dat$data, NULL, dat$compare, dat$index),
    "'model' must be a dust_generator")
  expect_error(
    particle_nofilter$new(dat$data, dat$model, dat$compare, TRUE),
    "'index' must be function if not NULL")
  expect_error(
    particle_nofilter$new(dat$data, dat$model, dat$compare, dat$index,
                          initial = 1),
    "'initial' must be function if not NULL")
  data_nested <- structure(dat$data,
                           class = c("particle_filter_data",
                                     "particle_filter_data_nested"))
  expect_error(
    particle_nofilter$new(data_nested, dat$model, dat$compare, dat$index),
    "nested mode not yet supported")
})


test_that("Can run mcmc with deterministic filter", {
  dat <- example_sir()
  control <- pmcmc_control(100, save_trajectories = TRUE)
  p <- particle_nofilter$new(dat$data, dat$model, dat$compare, dat$index)
  res <- pmcmc(dat$pars, p, NULL, control)
  expect_s3_class(res, "mcstate_pmcmc")
})


test_that("returning seed varies as model is run", {
  dat <- example_sir()
  control <- pmcmc_control(100, save_trajectories = TRUE)
  p <- particle_nofilter$new(dat$data, dat$model, dat$compare, dat$index,
                             seed = 1L)
  expect_equal(p$inputs()$seed, 1L)
  p$run()
  expect_type(p$inputs()$seed, "raw")
})


test_that("Can run deterministic filter without index", {
  dat <- example_sir()
  n_particles <- 42
  p1 <- particle_nofilter$new(dat$data, dat$model, dat$compare,
                              index = dat$index)

  compare2 <- function(state, ...) {
    dat$compare(state[5, , drop = FALSE], ...)
  }

  p2 <- particle_nofilter$new(dat$data, dat$model, compare2)

  set.seed(1)
  ll1 <- p1$run(save_history = TRUE)
  set.seed(1)
  ll2 <- p2$run(save_history = TRUE)
  expect_identical(ll1, ll2)

  expect_equal(dim(p1$history()), c(3, 1, 101))
  expect_equal(dim(p2$history()), c(5, 1, 101))
})
