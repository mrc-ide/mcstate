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
