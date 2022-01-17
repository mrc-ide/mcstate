context("pmcmc (parallel)")

test_that("basic parallel operation", {
  dat <- example_sir()
  n_particles <- 10
  n_steps <- 15
  n_chains <- 4

  filter <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                                index = dat$index, seed = 1L)

  ## TODO: I am not clear why 'use_parallel_seed' here was apparently
  ## optional; probably worth preserving?
  control_serial <- pmcmc_control(n_steps, n_chains = n_chains,
                                  progress = FALSE, use_parallel_seed = TRUE)
  cmp <- pmcmc(dat$pars, filter, control = control_serial)

  control_parallel <- pmcmc_control(n_steps, n_chains = n_chains,
                                    n_workers = 2L, n_threads_total = 2L,
                                    progress = FALSE, use_parallel_seed = TRUE)
  ans <- pmcmc(dat$pars, filter, control = control_parallel)

  expect_equal(cmp$pars, ans$pars)
  expect_equal(cmp, ans)
})


test_that("running pmcmc with progress = TRUE prints messages", {
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                           index = dat$index, seed = 1L)
  control <- pmcmc_control(10, n_chains = 2, n_workers = 2L, progress = TRUE)
  expect_message(
    pmcmc(dat$pars, p, control = control),
    "Finished 20 steps in ")
})


test_that("make seeds with integer returns length-32 raws", {
  algorithm <- "xoshiro256plus"
  s <- make_seeds(2, 1, algorithm)
  expect_length(s, 2L)
  expect_setequal(names(s[[1]]), c("dust", "r"))
  expect_length(s[[1]]$dust, 32L)
  expect_is(s[[1]]$dust, "raw")
  expect_identical(make_seeds(2, 1, algorithm), s)
  expect_identical(make_seeds(4, 1, algorithm)[1:2], s)
})


test_that("make seeds with long raw retains size", {
  algorithm <- "xoshiro256plus"
  seed <- dust::dust_rng$new(NULL, n_streams = 3)$state()
  s <- make_seeds(2L, seed, algorithm)
  expect_length(s, 2L)
  expect_setequal(names(s[[1]]), c("dust", "r"))
  expect_identical(s[[1]]$dust, seed)
  expect_length(s[[2]]$dust, length(seed))
})


test_that("construct parallel filter data", {
  expect_equal(
    pmcmc_parallel_progress_data(c("pending", "pending"), c(0, 0), 100),
    list(steps = 0,
         tokens = list(bar_overall = "  ", p_running = "")))
  expect_equal(
    pmcmc_parallel_progress_data(rep("pending", 4), rep(0, 4), 100),
    list(steps = 0,
         tokens = list(bar_overall = "    ", p_running = "")))
  expect_equal(
    pmcmc_parallel_progress_data(c("running", "pending", "pending", "pending"),
                                 c(50, 0, 0, 0), 100),
    list(steps = 50,
         tokens = list(bar_overall = "+   ", p_running = " 50%")))
  expect_equal(
    pmcmc_parallel_progress_data(
      c("done", "running", "running", "pending", "pending"),
      c(100, 12, 67, 0, 0),
      100),
    list(steps = 179,
         tokens = list(bar_overall = "#++  ", p_running = " 12%  67%")))
  expect_equal(
    pmcmc_parallel_progress_data(c("done", "done", "running"), c(100, 100, 95),
                                 100),
    list(steps = 295,
         tokens = list(bar_overall = "##+", p_running = " 95%")))
  expect_equal(
    pmcmc_parallel_progress_data(rep("done", 3), rep(100, 3), 100),
    list(steps = 300,
         tokens = list(bar_overall = "###", p_running = "")))
})


test_that("progress bar is a noop when progress = FALSE", {
  control <- pmcmc_control(100, n_workers = 2, n_chains = 3, progress = FALSE)
  p <- pmcmc_parallel_progress(control, force = TRUE)
  Sys.sleep(0.2)
  expect_silent(p(rep("pending", 3), rep(0, 3)))
  expect_silent(p(c("running", "pending", "pending"), c(100, 0, 0)))
  expect_silent(p(c("done", "done", "done"), c(100, 100, 100)))
})


test_that("progress bar creates progress_bar when progress = TRUE", {
  control <- pmcmc_control(100, n_workers = 2, n_chains = 2, progress = TRUE)
  status <- rep("pending", 2)
  steps <- rep(0, 2)
  p <- pmcmc_parallel_progress(control, status, steps, force = TRUE)

  Sys.sleep(0.2)
  expect_message(
    p(c("running", "pending"), c(50, 0)),
    "\\[.\\] \\[\\+ \\] ETA .* \\| 00:00:[0-9]{2} so far \\( 50%\\)")
  expect_message(
    p(c("running", "running"), c(90, 10)),
    "\\[.\\] \\[\\+\\+\\] ETA .* \\| 00:00:[0-9]{2} so far \\( 90%  10%\\)")
  expect_message(
    p(c("done", "running"), c(100, 33)),
    "\\[.\\] \\[\\#\\+\\] ETA .* \\| 00:00:[0-9]{2} so far \\( 33%\\)")
})


test_that("basic parallel operation nested", {
  dat <- example_sir_shared()
  n_particles <- 42
  n_steps <- 30
  n_chains <- 4

  proposal_fixed <- matrix(0.00026)
  proposal_varied <- matrix(0.00057)

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = proposal_fixed, proposal_varied = proposal_varied)

  filter <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                                index = dat$index, seed = 1L)
  control_parallel <- pmcmc_control(n_steps, n_chains = n_chains,
                                    n_workers = 2L)
  ans <- pmcmc(pars, filter, control = control_parallel)

  control_serial <- pmcmc_control(n_steps, n_chains = n_chains,
                                  progress = FALSE, use_parallel_seed = TRUE)
  cmp <- pmcmc(dat$pars, filter, control = control_serial)

  expect_equal(cmp$pars, ans$pars)
  expect_equal(cmp, ans)
})


test_that("pmcmc can save files in given place", {
  dat <- example_sir()
  n_particles <- 10
  n_steps <- 15
  n_chains <- 2
  filter <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                                index = dat$index, seed = 1L)

  ## TODO: I am not clear why 'use_parallel_seed' here was apparently
  path <- tempfile()
  control <- pmcmc_control(n_steps, n_chains = n_chains, n_workers = 2L,
                           n_threads_total = 2L,
                           progress = FALSE, use_parallel_seed = TRUE,
                           path = path)
  ans <- pmcmc(dat$pars, filter, control = control)
  expect_true(file.exists(path))
  expect_true(file.exists(file.path(path, "results_1.rds")))
  expect_true(file.exists(file.path(path, "inputs.rds")))
})
