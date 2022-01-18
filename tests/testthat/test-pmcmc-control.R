context("pmcmc (control)")


test_that("Don't allow fewer chains than workers", {
  expect_error(
    pmcmc_control(100, n_chains = 2, n_workers = 5, n_threads_total = 5),
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


test_that("specify save_restart", {
  expect_null(pmcmc_control(100, save_restart = NULL)$save_restart)
  expect_equal(pmcmc_control(100, save_restart = 1)$save_restart, 1)
  expect_equal(pmcmc_control(100, save_restart = c(10, 20))$save_restart,
               c(10, 20))
  expect_error(
    pmcmc_control(100, save_restart = c(20, 10)),
    "'save_restart' must be strictly increasing")
})

test_that("integer step ratio", {
  expect_error(pmcmc_control(1, nested_step_ratio = 1.2), "must be an integer")
  expect_error(pmcmc_control(1, nested_step_ratio = 1 / 1.2), "must be")
  expect_silent(pmcmc_control(1, nested_step_ratio = 3))
  expect_silent(pmcmc_control(1, nested_step_ratio = 1 / 3))
})


test_that("filter on generation - no filter", {
  dat <- pmcmc_filter_on_generation(100, NULL, NULL)
  expect_equal(dat, list(n_burnin = 0, n_steps_retain = 100, n_steps_every = 1))
  steps <- seq(dat$n_burnin + 1, by = dat$n_steps_every,
               length.out = dat$n_steps_retain)
  expect_equal(steps, 1:100)
  i <- seq_len(100)
  expect_equal(
    which(i >= dat$n_burnin & (i - dat$n_burnin - 1) %% dat$n_steps_every == 0),
    steps)

  expect_equal(
    dat$n_burnin + (dat$n_steps_retain - 1) * dat$n_steps_every + 1,
    100)
})


test_that("filter on generation - burnin and filter", {
  dat <- pmcmc_filter_on_generation(100, 40, 20)
  expect_equal(dat, list(n_burnin = 42, n_steps_retain = 20, n_steps_every = 3))
  steps <- seq(dat$n_burnin + 1, by = dat$n_steps_every,
               length.out = dat$n_steps_retain)
  expect_equal(steps, seq(43, 100, by = 3))
  i <- seq_len(100)
  expect_equal(
    which(i >= dat$n_burnin & (i - dat$n_burnin - 1) %% dat$n_steps_every == 0),
    steps)

  expect_equal(
    dat$n_burnin + (dat$n_steps_retain - 1) * dat$n_steps_every + 1,
    100)
})


test_that("prevent invalid burnin and filter", {
  expect_error(
    pmcmc_filter_on_generation(10, 100, 5),
    "'n_burnin' cannot be greater than or equal to 'n_steps'")
  expect_error(
    pmcmc_filter_on_generation(100, 100, 5),
    "'n_burnin' cannot be greater than or equal to 'n_steps'")
  expect_error(
    pmcmc_filter_on_generation(100, 10, 500),
    "'n_steps_retain' is too large, max possible is 90 but given 500")
  expect_error(
    pmcmc_filter_on_generation(100, 10, 75),
    "'n_steps_retain' is too large to skip any samples,")
})


test_that("control can detect corruption", {
  control <- pmcmc_control(100, n_steps_retain = 15, n_burnin = 5)
  control$n_steps <- 30
  expect_error(
    pmcmc_check_control(control),
    "Corrupt pmcmc_control (n_steps/n_steps_retain/n_burnin)",
    fixed = TRUE)
})


test_that("Informational notices in control", {
  expect_message(
    ctrl <- pmcmc_control(10, path = "location"),
    "'path' given without n_workers has no effect and is ignored")
})
