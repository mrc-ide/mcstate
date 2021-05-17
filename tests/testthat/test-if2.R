context("IF2")

test_that("Can run IF2", {
  dat <- example_sir()

  pars <- if2_parameters$new(
            list(if2_parameter("beta", 0.15, min = 0, max = 1),
                 if2_parameter("gamma", 0.05, min = 0, max = 1)))

  iterations <- 50
  cooling_target <- 0.5
  n_par_sets <- 20
  control <- if2_control(pars_sd = list("beta" = 0.02, "gamma" = 0.02),
                         iterations = iterations,
                         n_par_sets = n_par_sets,
                         cooling_target = cooling_target,
                         progress = FALSE)

  filter <- if2$new(pars, dat$data, dat$model, dat$compare, NULL,
                    dat$index, control)

  set.seed(1)
  filter$run()

  expect_equal(length(filter$log_likelihood()), iterations)
  expect_equal(dim(filter$pars_series()),
               c(length(pars$names()), n_par_sets, iterations))

  # Test plots can run
  filter$plot() # LL plot (default)
  filter$plot("beta")
  filter$plot("gamma")

  # LL from particle filter runs possible
  n_particles <- 50
  ll_samples <- filter$sample(n_particles)
  expect_equal(length(ll_samples), n_par_sets)
})

test_that("IF2 won't run with mismatched parameter names", {
  dat <- example_sir()

  pars <- if2_parameters$new(
            list(if2_parameter("b", 0.15, min = 0, max = 1),
                 if2_parameter("gamma", 0.05, min = 0, max = 1)))

  iterations <- 50
  cooling_target <- 0.5
  n_par_sets <- 20
  control <- if2_control(pars_sd = list("beta" = 0.02, "gamma" = 0.02),
                         iterations = iterations,
                         n_par_sets = n_par_sets,
                         cooling_target = cooling_target,
                         progress = FALSE)

  expect_error(if2$new(pars, dat$data, dat$model, dat$compare, NULL,
                       dat$index, control),
               "'{b}' must be in control$pars_sd",
               fixed = TRUE)
})

test_that("IF2 inputs must be of the correct type", {
  dat <- example_sir()

  pars <- if2_parameters$new(
            list(if2_parameter("beta", 0.15, min = 0, max = 1),
                 if2_parameter("gamma", 0.05, min = 0, max = 1)))

  iterations <- 50
  cooling_target <- 0.5
  n_par_sets <- 20
  control <- if2_control(pars_sd = list("beta" = 0.02, "gamma" = 0.02),
                         iterations = iterations,
                         n_par_sets = n_par_sets,
                         cooling_target = cooling_target,
                         progress = FALSE)

  expect_error(if2$new(pars, dat$data, "dat$model", dat$compare, NULL,
                       dat$index, control),
                       "'model' must be a dust_generator")
  expect_error(if2$new(pars, dat$data, dat$model, dat$compare, NULL,
                       c(1, 2), control),
                       "'index' must be function if not NULL")
  expect_error(if2$new(pars, dat$data, dat$model, "dat$compare", NULL,
                        dat$index, control),
                       "'compare' must be a function")
})

test_that("Can't get IF2 results before object has been run", {
  dat <- example_sir()

  pars <- if2_parameters$new(
            list(if2_parameter("beta", 0.15, min = 0, max = 1),
                 if2_parameter("gamma", 0.05, min = 0, max = 1)))

  iterations <- 50
  cooling_target <- 0.5
  n_par_sets <- 20
  control <- if2_control(pars_sd = list("beta" = 0.02, "gamma" = 0.02),
                         iterations = iterations,
                         n_par_sets = n_par_sets,
                         cooling_target = cooling_target,
                         progress = FALSE)

  filter <- if2$new(pars, dat$data, dat$model, dat$compare, NULL,
                    dat$index, control)
  expect_error(filter$log_likelihood(), "IF2 must be run first")
  expect_error(filter$pars_series(), "IF2 must be run first")
  expect_error(filter$plot(), "IF2 must be run first")
  expect_error(filter$sample(100L), "IF2 must be run first")
})

test_that("stop inference when likelihood is impossible", {
  dat <- example_sir()
  steps <- nrow(dat$data) + 1

  compare <- function(state, observed, pars) {
    ret <- dat$compare(state, observed, pars)
    if (observed$incidence > 15) {
      ret[] <- -Inf
    }
    ret
  }

  pars <- if2_parameters$new(
            list(if2_parameter("beta", 0.15, min = 0, max = 1),
                 if2_parameter("gamma", 0.05, min = 0, max = 1)))

  iterations <- 50
  cooling_target <- 0.5
  n_par_sets <- 20
  control <- if2_control(pars_sd = list("beta" = 0.02, "gamma" = 0.02),
                         iterations = iterations,
                         n_par_sets = n_par_sets,
                         cooling_target = cooling_target,
                         progress = FALSE)

  filter <- if2$new(pars, dat$data, dat$model, compare, NULL,
                    dat$index, control)
  filter$run()

  expect_equal(filter$log_likelihood(), rep(-Inf, iterations))
})
