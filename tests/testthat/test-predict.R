context("predict")

test_that("can run a prediction from a mcmc run", {
  dat <- example_sir()

  n_particles <- 100
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)

  set.seed(1)
  control <- pmcmc_control(30, save_state = TRUE, save_trajectories = TRUE)
  results <- pmcmc(dat$pars, p, control = control)

  steps <- seq(tail(dat$data$time_end, 1), by = 4, length.out = 26)
  y <- pmcmc_predict(results, steps)

  expect_equal(y$time, steps)
  expect_equal(y$rate, 4)
  expect_equal(y$predicted, rep(TRUE, length(steps)))

  expect_equal(dim(y$state), c(3, 30, length(steps)))

  ## Start from correct point
  expect_equal(y$state[, , 1], results$state[1:3, ])

  ## Check predictions are reasonable:
  expect_true(all(diff(t(y$state[1, , ])) <= 0))
  expect_true(all(diff(t(y$state[3, , ])) >= 0))
})


test_that("Can rehydrate mcmc results and predict", {
  results <- example_sir_pmcmc()$pmcmc

  steps <- seq(results$predict$time, by = 4, length.out = 26)
  y1 <- pmcmc_predict(results, steps, seed = 1L)

  results <- unserialize(serialize(results, NULL))
  expect_identical(pmcmc_predict(results, steps, seed = 1L), y1)
})


test_that("Can combine runs with predictions", {
  results <- example_sir_pmcmc()$pmcmc

  steps <- seq(results$predict$time, by = 4, length.out = 26)
  y1 <- pmcmc_predict(results, steps, prepend_trajectories = TRUE, seed = 1L)
  y2 <- pmcmc_predict(results, steps, seed = 1L)

  expect_is(y1, "mcstate_trajectories")
  expect_equal(y1$rate, 4)
  expect_equal(y1$time, seq(0, 500, by = 4))
  expect_equal(y1$predicted, rep(c(FALSE, TRUE), c(101, 25)))
  expect_identical(y1$state[, , 101:126], y2$state)
  expect_equal(y1$state[, , 1:101], results$trajectories$state)
})


test_that("Do not run a predict with no state", {
  dat <- example_uniform()
  control <- pmcmc_control(1000, save_state = FALSE)
  res <- pmcmc(dat$pars, dat$filter, control = control)
  expect_error(
    pmcmc_predict(res, 0:10),
    "mcmc was run with return_state = FALSE, can't predict")
})


test_that("Require at least two times", {
  results <- example_sir_pmcmc()$pmcmc
  expect_error(
    pmcmc_predict(results, integer(0)),
    "At least two times required for predict")
  expect_error(
    pmcmc_predict(results, results$predict$time),
    "At least two times required for predict")
})


test_that("Require start times to be correct", {
  results <- example_sir_pmcmc()$pmcmc
  steps <- seq(results$predict$time, by = 4, length.out = 26)
  expect_error(
    pmcmc_predict(results, steps + 1),
    "Expected times[1] to be 400",
    fixed = TRUE)
  expect_error(
    pmcmc_predict(results, steps - 1),
    "Expected times[1] to be 400",
    fixed = TRUE)
})


test_that("S3 method works as expected", {
  results <- example_sir_pmcmc()$pmcmc
  steps <- seq(results$predict$time, by = 4, length.out = 26)
  set.seed(1)
  y1 <- predict(results, steps)
  set.seed(1)
  y2 <- pmcmc_predict(results, steps)
  expect_identical(y1, y2)
})


test_that("seed = NULL respects R's seed", {
  results <- example_sir_pmcmc()$pmcmc
  steps <- seq(results$predict$time, by = 4, length.out = 26)
  set.seed(1)
  y1 <- pmcmc_predict(results, steps)
  y2 <- pmcmc_predict(results, steps)
  set.seed(1)
  y3 <- pmcmc_predict(results, steps)
  expect_false(identical(y1, y2))
  expect_identical(y1, y3)
})


test_that("can't prepend history if trajectories not saved", {
  results <- example_sir_pmcmc()$pmcmc
  results$trajectories <- NULL

  steps <- seq(results$predict$time, by = 4, length.out = 26)
  expect_error(
    pmcmc_predict(results, steps, prepend_trajectories = TRUE),
    "mcmc was run with return_trajectories = FALSE, can't prepend trajectories")
})


test_that("names are copied from index into predictions", {
  results <- example_sir_pmcmc()$pmcmc
  steps <- seq(results$predict$time, by = 4, length.out = 26)
  y1 <- pmcmc_predict(results, steps, seed = 1L)

  names(results$predict$index) <- c("a", "b", "c")
  y2 <- pmcmc_predict(results, steps, seed = 1L)
  y3 <- pmcmc_predict(results, steps, prepend_trajectories = TRUE, seed = 1L)

  expect_null(rownames(y1$state))
  expect_equal(rownames(y2$state), c("a", "b", "c"))
  expect_equal(rownames(y3$state), c("a", "b", "c"))
})

test_that("can run a prediction from a nested mcmc run", {
  dat <- example_sir_shared()

  n_particles <- 100
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)

  set.seed(1)
  control <- pmcmc_control(30, save_state = TRUE, save_trajectories = TRUE)
  results <- pmcmc(dat$pars, p, control = control)

  steps <- seq(tail(dat$data$time_end, 1), by = 4, length.out = 26)
  y <- pmcmc_predict(results, steps)

  expect_equal(y$time, steps)
  expect_equal(y$rate, 4)
  expect_equal(y$predicted, rep(TRUE, length(steps)))

  expect_equal(dim(y$state), c(3, 2, 30, length(steps)))

  ## Check predictions are reasonable:
  expect_true(all(diff(t(y$state[1, 1, , ])) <= 0))
  expect_true(all(diff(t(y$state[3, 2, , ])) >= 0))

  expect_silent(pmcmc_predict(results, steps, prepend_trajectories = TRUE))
})


test_that("Can't predict from a continuous pmcmc", {
  dat <- example_continuous()
  n_particles <- 20
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L,
                           stochastic_schedule = dat$stochastic_schedule)
  control <- pmcmc_control(10, save_trajectories = TRUE, save_state = TRUE)
  results <- pmcmc(dat$pars, p, control = control)
  expect_error(
    pmcmc_predict(results, seq(1800, length.out = 10, by = 30)),
    "predict not (yet) possible with continuous models (mrc-3453)",
    fixed = TRUE)
})
