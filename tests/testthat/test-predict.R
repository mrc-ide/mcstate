context("predictn")

test_that("can run a prediction from a mcmc run", {
  dat <- example_sir()

  n_particles <- 100
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)

  set.seed(1)
  results <- pmcmc(dat$pars, p, 30, TRUE, TRUE)

  steps <- seq(tail(dat$data$step_end, 1), by = 4, length.out = 26)
  y <- pmcmc_predict(results, steps)

  expect_equal(dim(y), c(3, 31, length(steps)))

  ## Start from correct point
  expect_equal(y[, , 1], results$state[1:3, ])

  ## Check predictions are reasonable:
  expect_true(all(diff(t(y[1, , ])) <= 0))
  expect_true(all(diff(t(y[3, , ])) >= 0))
})
