context("pmcmc (control)")


test_that("Don't allow fewer chains than workers", {
  expect_error(
    pmcmc_control(100, n_chains = 2, n_workers = 5),
    "'n_chains' (2) is less than 'n_workers' (5)",
    fixed = TRUE)
})


test_that("Don't run with invalid n_threads", {
  expect_error(
    pmcmc_control(100, n_chains = 5, n_workers = 5, n_threads_total = 3),
    "'n_threads_total' (3) is less than 'n_workers' (5)",
    fixed = TRUE)
  expect_error(
    pmcmc_control(100, n_chains = 5, n_workers = 5, n_threads_total = 8),
    "'n_threads_total' (8) is not a multiple of 'n_workers' (5)",
    fixed = TRUE)
})


test_that("If not using workers, n_run_each has no effect", {
  expect_equal(
    pmcmc_control(100, n_workers = 1, n_steps_each = 10)$n_steps_each,
    100)
})


test_that("If not given, n_run_each is 10% of n_steps", {
  expect_equal(
    pmcmc_control(100, n_chains = 2, n_workers = 2)$n_steps_each,
    10)
  expect_equal(
    pmcmc_control(105, n_chains = 2, n_workers = 2)$n_steps_each,
    11)
})


test_that("specify save_restart", {
  expect_null(pmcmc_control(100, save_restart = NULL)$save_restart)
  expect_equal(pmcmc_control(100, save_restart = 1)$save_restart, 1)
  expect_equal(pmcmc_control(100, save_restart = c(10, 20))$save_restart,
               c(10, 20))
  expect_error(
    pmcmc_control(100, save_restart = c(20, 10)),
    "'save_restart' must be strictly increasing")
})
