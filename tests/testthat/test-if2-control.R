context("IF2 (control)")


test_that("Don't run with invalid pars_sd", {
  expect_error(
    if2_control(pars_sd = c("beta" = 0.02, "gamma" = 0.02),
                iterations = 100, n_par_sets = 100,
                cooling_target = 0.5, progress = TRUE),
    "'pars_sd' must be a list",
    fixed = TRUE)
  expect_error(
    if2_control(pars_sd = list("beta" = "1", "gamma" = 0.02),
                iterations = 100, n_par_sets = 100,
                cooling_target = 0.5, progress = TRUE),
    "Elements of 'pars_sd' must be in '{numeric}'",
    fixed = TRUE)
  expect_error(
    if2_control(pars_sd = list(0.02, 0.02),
                iterations = 100, n_par_sets = 100,
                cooling_target = 0.5, progress = TRUE),
    "'pars_sd' must be named",
    fixed = TRUE)
})

test_that("Don't run with invalid IF2 parameters", {
  expect_error(
    if2_control(pars_sd = list("beta" = 0.02, "gamma" = 0.02),
                iterations = -1, n_par_sets = 100,
                cooling_target = 0.5, progress = TRUE),
    "'iterations' must be at least 1",
    fixed = TRUE)
  expect_error(
    if2_control(pars_sd = list("beta" = 0.02, "gamma" = 0.02),
                iterations = 100, n_par_sets = -1,
                cooling_target = 0.5, progress = TRUE),
    "'n_par_sets' must be at least 1",
    fixed = TRUE)
  expect_error(
    if2_control(pars_sd = list("beta" = 0.02, "gamma" = 0.02),
                iterations = 100, n_par_sets = 100,
                cooling_target = 0, progress = TRUE),
    "'cooling_target' must be between 0 and 1 (non-inclusive)",
    fixed = TRUE)
  expect_error(
    if2_control(pars_sd = list("beta" = 0.02, "gamma" = 0.02),
                iterations = 100, n_par_sets = 100,
                cooling_target = 1, progress = TRUE),
    "'cooling_target' must be between 0 and 1 (non-inclusive)",
    fixed = TRUE)
})
