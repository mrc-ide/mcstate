context("forecast")

test_that("Sampling and forecasting from a grid search", {
  skip("redo")
  range <- data.frame(name = c("beta", "gamma"),
                      min = c(0.13, 0.05),
                      max = c(0.25, 0.15),
                      n = c(6, 9),
                      target = "pars_model",
                      stringsAsFactors = FALSE)

  dat <- example_sir()
  n_particles <- 100
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  state <- dat$y0
  n_sample_pars <- 10
  forecast_steps <- 0

  set.seed(1)
  grid_res <- grid_search(range, p)
  forecast_res <- forecast(grid_res,
                           filter = p,
                           n_sample_pars = n_sample_pars,
                           forecast_steps = forecast_steps)

  # check structure is as expected
  expect_is(forecast_res, "mcstate_forecast")
  expect_equal(dim(forecast_res$parameters), c(n_sample_pars, nrow(range)))
  expect_equal(names(forecast_res$parameters), range$name)
  # Number of parameter samples as expected
  expect_length(forecast_res$trajectories, n_sample_pars)
  for (i in seq_len(n_sample_pars)) {
    # 3 states, 101 steps
    expect_equal(dim(forecast_res$trajectories[[i]]), c(3, n_particles, 101))
  }

  # check that the input data is recovered
  # averages over all parameter samples and all particles
  grid <- expand.grid(x = seq_len(n_sample_pars), y = seq_len(n_particles))

  f <- function(i, j, k) {
    fit <- lm(forecast_res$trajectories[[i]][k, j, ] ~ dat$history[k, 1, ])
    summary(fit)$adj.r.squared
  }

  s_r2 <- mapply(f, grid$x, grid$y, 1)
  expect_gt(mean(s_r2), 0.99)

  # I is much more noisy - and a better test case
  i_r2 <- mapply(f, grid$x, grid$y, 2)
  expect_gt(mean(i_r2), 0.95)

  r_r2 <- mapply(f, grid$x, grid$y, 3)
  expect_gt(mean(r_r2), 0.99)

  plot(forecast_res, what = 2, data = dat$history[2, 1, ],
       title = "I", ylab = "I")

  # check that forecasting is possible
  forecast_steps <- 5
  forecast_res <- forecast(grid_res,
                           filter = p,
                           n_sample_pars = n_sample_pars,
                           forecast_steps = forecast_steps)
  for (i in seq_len(n_sample_pars)) {
    # 5 quantities, 101 steps
    expect_equal(dim(forecast_res$trajectories[[i]]), c(3, n_particles, 106))
  }

  plot(forecast_res, what = 3, data = dat$history[3, 1, ],
       title = "R", ylab = "R")

})

test_that("Sampling and forecasting from an MCMC", {
  skip("redo")
  range <- data.frame(name = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  lprior <- list("beta" = function(pars) log(1e-10),
                 "gamma" = function(pars) log(1e-10))
  proposal_kernel <- diag(nrow(range)) * 0.01^2
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- range$name

  dat <- example_sir()
  n_particles <- 20
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  n_mcmc <- 100
  n_chains <- 2
  n_sample_pars <- 10
  forecast_steps <- 0

  set.seed(1)
  mcmc_res <- pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                    n_chains = n_chains)
  forecast_res <- forecast(mcmc_res,
                           filter = p,
                           n_sample_pars = n_sample_pars,
                           forecast_steps = forecast_steps,
                           burn_in = 10)

  # check structure is as expected
  expect_is(forecast_res, "mcstate_forecast")
  expect_equal(dim(forecast_res$parameters), c(n_sample_pars, nrow(range)))
  expect_equal(names(forecast_res$parameters), range$name)
  # Number of parameter samples as expected
  expect_length(forecast_res$trajectories, n_sample_pars)
  for (i in seq_len(n_sample_pars)) {
    # 3 states, 101 steps
    expect_equal(dim(forecast_res$trajectories[[i]]), c(3, n_particles, 101))
  }

  # check that the input data is recovered
  # averages over all parameter samples and all particles
  grid <- expand.grid(x = seq_len(n_sample_pars), y = seq_len(n_particles))

  f <- function(i, j, k) {
    fit <- lm(forecast_res$trajectories[[i]][k, j, ] ~ dat$history[k, 1, ])
    summary(fit)$adj.r.squared
  }

  s_r2 <- mapply(f, grid$x, grid$y, 1)
  expect_gt(mean(s_r2), 0.99)

  # I is much more noisy - and a better test case
  i_r2 <- mapply(f, grid$x, grid$y, 2)
  expect_gt(mean(i_r2), 0.95)

  r_r2 <- mapply(f, grid$x, grid$y, 3)
  expect_gt(mean(r_r2), 0.99)

  plot(forecast_res, what = 2, data = dat$history[2, 1, ],
       title = "I", ylab = "I")

  # check that forecasting is possible
  forecast_steps <- 5
  forecast_res <- forecast(mcmc_res,
                           filter = p,
                           n_sample_pars = n_sample_pars,
                           forecast_steps = forecast_steps)
  for (i in seq_len(n_sample_pars)) {
    # 5 quantities, 101 steps
    expect_equal(dim(forecast_res$trajectories[[i]]), c(3, n_particles, 106))
  }

  plot(forecast_res, what = 3, data = dat$history[2, 1, ],
       title = "R", ylab = "R")

})


test_that("Can't plot outside of range", {
  x <- list(trajectories = list(matrix(NA, 10)))
  expect_error(
    plot.mcstate_forecast(x, what = 0),
    "'what' must be a valid index for a partition")
  expect_error(
    plot.mcstate_forecast(x, what = 11),
    "'what' must be a valid index for a partition")
})
