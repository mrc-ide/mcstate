context("grid_search")

test_that("Simple grid search with SIR model", {

  range <- data_frame(name = c("beta", "gamma"),
                      min = c(0.15, 0.05),
                      max = c(0.3, 0.15),
                      n = c(9, 9),
                      target = "model_data")

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
  range <- data_frame(name = c("beta", "gamma"),
                      min = c(0.1, 0),
                      max = c(0.3, 0.2),
                      n = c(3, 3),
                      target = "model_data")

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

test_that("Start date can be sampled", {
  range <- data.frame(name = c("beta", "step_start"),
                      min = c(0.1, 0),
                      max = c(0.3, 100),
                      n = c(3, 3),
                      target = c("model_data", "step_start"),
                      stringsAsFactors = FALSE)

  dat <- example_sir()
  data <- dat$data
  offset <- 100
  data[c("step_start", "step_end")] <-
    data[c("step_start", "step_end")] + offset
  data$step_start[[1]] <- 0

  p <- particle_filter$new(data, dat$model, dat$compare)
  n_particles <- 100

  set.seed(1)
  grid_res <- grid_search(range, p, n_particles)
  expect_true(grid_res$renorm_mat_ll[2, 3] == max(grid_res$renorm_mat_ll))
})

test_that("pars_compare can be sampled", {
  range <- data.frame(name = c("beta", "exp_noise"),
                      min = c(0.1, 1e3),
                      max = c(0.3, 1e6),
                      n = c(3, 3),
                      target = c("model_data", "pars_compare"),
                      stringsAsFactors = FALSE)

  dat <- example_sir()

  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  n_particles <- 100

  set.seed(1)
  grid_res <- grid_search(range, p, n_particles)
  expect_true(grid_res$renorm_mat_ll[2, 1] == max(grid_res$renorm_mat_ll))
})

test_that("grid_search_validate_range - happy path", {
  range <- data_frame(name = c("a", "b"),
                      min = c(1, 0),
                      max = c(2, 10),
                      n = c(4, 5),
                      target = "model_data")
  res <- grid_search_validate_range(range)

  expect_setequal(
    names(res),
    c("range", "variables", "expanded", "index"))
  expect_equal(
    res$range, range)
  expect_equal(
    res$variables,
    list(a = seq(1, 2, length.out = 4),
         b = seq(0, 10, length.out = 5)))
  expect_equal(
    res$expanded,
    expand.grid(a = res$variables$a, b = res$variables$b))
  expect_equal(
    res$index,
    list(step_start = integer(0),
         model_data = 1:2,
         pars_compare = integer(0)))
})


test_that("grid_search_validate_range requires columns", {
  expect_error(
    grid_search_validate_range(
      data_frame(name = c("a", "b"), min = 1, max = 2, n = 2)),
    "Missing columns from 'range': 'target'")
  expect_error(
    grid_search_validate_range(
      data_frame(name = c("a", "b"), n = 2)),
    "Missing columns from 'range': 'min', 'max', 'target'")
})


test_that("grid_search_validate_range validates target", {
  expect_error(
    grid_search_validate_range(
      data_frame(name = c("a", "b"), min = 1, max = 2, n = 2,
                 target = "somewhere")),
    paste("Invalid target 'somewhere': must be one of 'step_start',",
          "'model_data', 'pars_compare'"))
  expect_error(
    grid_search_validate_range(
      data_frame(name = c("a", "b"), min = 1, max = 2, n = 2,
                 target = c("x", "y"))),
    paste("Invalid target 'x', 'y': must be one of 'step_start',",
          "'model_data', 'pars_compare'"))
})


test_that("grid_search_validate_range allows at most one step_start target", {
  expect_error(
    grid_search_validate_range(
      data_frame(name = c("a", "b"), min = 1, max = 2, n = 2,
                 target = "step_start")),
    "At most one target may be 'step_start'")
})


test_that("grid_search_validate_range requires two variables", {
  range <- data_frame(name = c("a", "b", "c"),
                      min = 1, max = 2, n = 10, target = "model_data")
  expect_error(
    grid_search_validate_range(range[1, ]),
    "Expected exactly two rows in 'range'")
  expect_error(
    grid_search_validate_range(range),
    "Expected exactly two rows in 'range'")
})


test_that("grid_search_validate_range requires unique names", {
  range <- data_frame(name = c("a", "a"),
                      min = 1, max = 2, n = 10, target = "model_data")
  expect_error(
    grid_search_validate_range(range),
    "Duplicate 'name' entries not allowed in 'range'")
})
