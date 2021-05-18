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

  filter <- particle_filter$new(dat$data, dat$model, 1L, dat$compare,
                                dat$index)

  obj <- if2(pars, filter, control)

  expect_s3_class(obj, "if2_fit")
  expect_length(obj$result$log_likelihood, iterations)
  expect_equal(dim(obj$result$pars),
               c(2, n_par_sets, iterations))

  ## LL from particle filter runs possible
  n_particles <- 50
  ll_samples <- if2_sample(obj, n_particles)
  expect_equal(length(ll_samples), n_par_sets)
})

test_that("IF2 converges on the correct ll for the volatility model", {
  dat <- example_volatility()

  # Start parameters far away from correct values
  set.seed(1)
  transform <- function(x) {
    c(as.list(x), list(compare = list(gamma = 1, tau = 1)))
  }
  pars <- if2_parameters$new(
    list(if2_parameter("alpha", 5, min = 0, max = Inf),
         if2_parameter("sigma", 5, min = 0, max = Inf)),
    transform)

  iterations <- 50
  cooling_target <- 0.5
  n_par_sets <- 100
  control <- if2_control(pars_sd = list("alpha" = 0.02, "sigma" = 0.02),
                         iterations = iterations,
                         n_par_sets = n_par_sets,
                         cooling_target = cooling_target,
                         progress = FALSE)

  filter <- particle_filter$new(dat$data, dat$model, 1L, dat$compare)
  obj <- if2(pars, filter, control)

  filter_ll <- if2_sample(obj, 100)
  expect_equal(mean(filter_ll),
               dat$kalman_filter(dat$pars, dat$data),
               sd(filter_ll))
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

  filter <- particle_filter$new(dat$data, dat$model, 1L, dat$compare,
                                dat$index)

  expect_error(if2(pars, filter, control),
               "'{b}' must be in control$pars_sd",
               fixed = TRUE)
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
  n_par_sets <- 5
  control <- if2_control(pars_sd = list("beta" = 0.02, "gamma" = 0.02),
                         iterations = iterations,
                         n_par_sets = n_par_sets,
                         cooling_target = cooling_target,
                         progress = FALSE)

  filter <- particle_filter$new(dat$data, dat$model, 1L, compare,
                                dat$index)

  res <- if2(pars, filter, control)

  expect_equal(res$result$log_likelihood, rep(-Inf, iterations))
})
