context("grid_search")

test_that("Simple grid search with SIR model", {

  range <- data.frame(name = c("beta", "gamma"),
                      min = c(0.15, 0.05),
                      max = c(0.3, 0.15),
                      n = c(9, 9),
                      target = "model_data",
                      stringsAsFactors = FALSE)

  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  n_particles <- 100

  res <- grid_search(range, p, n_particles)

  expect_is(res, "mcstate_scan")

  beta_grid <- seq(range$min[1], range$max[1], length.out = range$n[1])
  gamma_grid <- seq(range$min[2], range$max[2], length.out = range$n[2])
  expect_equal(res$x$beta, beta_grid)
  expect_equal(res$y$gamma, gamma_grid)
  expect_equal(dim(res$renorm_mat_ll), dim(res$mat_log_ll))
  expect_equal(dim(res$renorm_mat_ll), c(length(beta_grid), length(gamma_grid)))
  expect_true(all(res$renorm_mat_ll <= 1 & res$renorm_mat_ll >= 0))

  plot(res)
  plot(res, what = "probability")
})

test_that("SIR model parameters are can be inferred correctly", {
  range <- data.frame(name = c("beta", "gamma"),
                      min = c(0.1, 0),
                      max = c(0.3, 0.2),
                      n = c(3, 3),
                      target = "model_data",
                      stringsAsFactors = FALSE)

  # NB: in the example model, default
  # * beta = 0.2 (transmission)
  # * g = 0.1 (recovery)
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  n_particles <- 100

  res <- grid_search(range, p, n_particles)

  # Correct parameter estimate is the highest
  expect_true(res$renorm_mat_ll[2, 2] == max(res$renorm_mat_ll))

})
