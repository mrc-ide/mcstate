context("grid_search")

test_that("Simple grid search with SIR model", {
  
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
  
  res <- grid_search(state, range, p, n_particles)
  
  expect_is(res, "grid_scan")
  #expect_true("inputs" %in% names(scan_results))
  #expect_setequal(names(scan_results$inputs),
  #                c("model", "model_params", "pars_obs", "data"))
  
  beta_grid <- seq(range$min[1], range$max[1], length.out = range$n[1])
  gamma_grid <- seq(range$min[2], range$max[2], length.out = range$n[2])
  expect_equal(res$x$beta, beta_grid)
  expect_equal(res$y$gamma, gamma_grid)
  expect_equal(dim(res$renorm_mat_LL), dim(res$mat_log_ll))
  expect_equal(dim(res$renorm_mat_LL), c(length(beta_grid), length(gamma_grid)))
  expect_true(all(res$renorm_mat_LL <= 1 & res$renorm_mat_LL >= 0))
  
  plot(res)
  plot(res, what='probability')
})

test_that("SIR model parameters are can be inferred correctly", {
  range <- data.frame(name = c("beta", "gamma"),
                      min = c(0.1, 0),
                      max = c(0.3, 0.2),
                      n = c(3, 3),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  
  # NB: in the example model, beta = 0.2 (transmission), g = 0.1 (recovery)
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  state <- dat$y0
  n_particles <- 100
  
  res <- grid_search(state, range, p, n_particles)
  
  # Correct parameter estimate is the highest
  expect_true(res$renorm_mat_LL[2,2] == max(res$renorm_mat_LL))
  
})
