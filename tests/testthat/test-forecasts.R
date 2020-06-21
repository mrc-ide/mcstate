context("Sampling and forecasts")

test_that("Sampling and forecasting from a grid search", {
  range <- data.frame(name = c("beta", "gamma"),
                      min = c(0.13, 0.05),
                      max = c(0.25, 0.15),
                      n = c(6, 9),
                      target = "model_data",
                      stringsAsFactors = FALSE)

  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  state <- dat$y0
  n_particles <- 100
  n_sample_pairs <- 10
  forecast_steps <- 0

  set.seed(1)
  grid_res <- grid_search(range, p, n_particles)
  forecast_res <- sample_grid_scan(grid_res, p, n_sample_pairs, n_particles,
                                   forecast_steps)

  # check structure is as expected
  expect_is(forecast_res, "mcstate_forecast")
  expect_equal(dim(forecast_res$parameters), c(n_sample_pairs, nrow(range)))
  expect_equal(names(forecast_res$parameters), range$name)
  # Number of parameter samples as expected
  expect_length(forecast_res$trajectories, n_sample_pairs)
  for (i in seq_len(n_sample_pairs)) {
    # 3 states, 101 steps
    expect_equal(dim(forecast_res$trajectories[[i]]), c(3, n_particles, 101))
  }

  # check that the input data is recovered
  # averages over all parameter samples and all particles
  grid <- expand.grid(x = seq_len(n_sample_pairs), y = seq_len(n_particles))
  s_r2 <- purrr::map2_dbl(.x = grid$x, .y = grid$y, .f = function(i, j) {
    fit <- lm(forecast_res$trajectories[[i]][1, j, ] ~ dat$history[1, 1, ])
    summary(fit)$adj.r.squared
  })
  expect_gt(mean(s_r2), 0.99)
  i_r2 <- purrr::map2_dbl(.x = grid$x, .y = grid$y, .f = function(i, j) {
    fit <- lm(forecast_res$trajectories[[i]][2, j, ] ~ dat$history[2, 1, ])
    summary(fit)$adj.r.squared
  })
  expect_gt(mean(i_r2), 0.95) # I is much more noisy - and a better test case
  r_r2 <- purrr::map2_dbl(.x = grid$x, .y = grid$y, .f = function(i, j) {
    fit <- lm(forecast_res$trajectories[[i]][3, j, ] ~ dat$history[3, 1, ])
    summary(fit)$adj.r.squared
  })
  expect_gt(mean(r_r2), 0.99)

  plot(dat$history[2, 1, ], type = "l", xlab = "day", ylab = "I", lwd = 2)
  for (j in seq_len(n_particles)) {
    lines(forecast_res$trajectories[[1]][2, j, ], type = "l", lty = 2,
          col = rgb(139, 0, 0, maxColorValue = 255, alpha = 30))
  }

  # check that forecasting is possible
  forecast_steps <- 5
  forecast_res <- sample_grid_scan(grid_res, p, n_sample_pairs, n_particles,
                                   forecast_steps)
  for (i in seq_len(n_sample_pairs)) {
    # 5 quantities, 101 steps
    expect_equal(dim(forecast_res$trajectories[[i]]), c(3, n_particles, 106))
  }

})
