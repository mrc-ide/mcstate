context("Sampling and forecasts")

test_that("Sampling and forecasting from a grid search", {
  range <- data.frame(name = c("beta", "gamma"),
                      min = c(0.13, 0.05),
                      max = c(0.25, 0.15),
                      n = c(6, 9),
                      target = "pars_model",
                      stringsAsFactors = FALSE)

  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  state <- dat$y0
  n_particles <- 100
  n_sample_pairs <- 10
  forecast_steps <- 0

  grid_res <- grid_search(state, range, p, n_particles)
  forecast_res <- sample_grid_scan(grid_res, p, n_sample_pairs, n_particles, forecast_steps)

  # check structure is as expected
  expect_is(forecast_res, "mcstate_forecast")
  expect_equal(dim(forecast_res$parameters), c(n_sample_pairs, nrow(range)))
  expect_equal(names(forecast_res$parameters), range$name)
  expect_length(forecast_res$trajectories, n_sample_pairs) # Number of parameter samples as expected
  for (i in 1:n_sample_pairs) {
    expect_equal(dim(forecast_res$trajectories[[i]]), c(5, n_particles, 101)) # 5 quantities, 101 steps
  }
  # check that the input data is recovered
  # averages over all parameter samples and all particles
  # TODO: would a difference be better
  grid <- expand.grid(x=1:n_sample_pairs, y=1:n_particles)
  S_r2 <- purrr::map2_dbl(.x = grid$x, .y=grid$y, .f = function(i,j) {
                    summary(lm(forecast_res$trajectories[[i]][1,j,] ~ dat$y[,'S']))$adj.r.squared})
  expect_gt(mean(S_r2), 0.99)
  I_r2 <- purrr::map2_dbl(.x = grid$x, .y=grid$y, .f = function(i,j) {
    summary(lm(forecast_res$trajectories[[i]][2,j,] ~ dat$y[,'I']))$adj.r.squared})
  expect_gt(mean(I_r2), 0.95) # I is much more noisy as it has a lower value throughout
  R_r2 <- purrr::map2_dbl(.x = grid$x, .y=grid$y, .f = function(i,j) {
    summary(lm(forecast_res$trajectories[[i]][3,j,] ~ dat$y[,'R']))$adj.r.squared})
  expect_gt(mean(R_r2), 0.99)

  # check that forecasting is possible
  forecast_steps <- 5
  forecast_res <- sample_grid_scan(grid_res, p, n_sample_pairs, n_particles, forecast_steps)
  for (i in 1:n_sample_pairs) {
    expect_equal(dim(forecast_res$trajectories[[i]]), c(5, n_particles, 106)) # 5 quantities, 101 steps
  }

})