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
  control_serial <- pmcmc_control(n_steps, n_chains = n_chains, progress = FALSE,
                                  use_parallel_seed = TRUE)
  cmp <- pmcmc(dat$pars, filter, control = control_serial)

  control_parallel <- pmcmc_control(n_steps, n_chains = n_chains, n_workers = 2L,
                                    n_threads_total = 2L,
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
  skip("rewrite")
  expect_equal(
    pmcmc_parallel_progress_data(list(NULL, NULL), 100),
    list(n = 0,
         tokens = list(bar_overall = "  ", p_running = ""),
         result = FALSE))
  expect_equal(
    pmcmc_parallel_progress_data(list(NULL, NULL, NULL, NULL), 100),
    list(n = 0,
         tokens = list(bar_overall = "    ", p_running = ""),
         result = FALSE))
  expect_equal(
    pmcmc_parallel_progress_data(list(list(step = 50, finished = FALSE),
                                      NULL, NULL, NULL), 100),
    list(n = 50,
         tokens = list(bar_overall = "+   ", p_running = "50%"),
         result = FALSE))
  expect_equal(
    pmcmc_parallel_progress_data(list(list(step = 100, finished = TRUE),
                                      list(step = 12, finished = FALSE),
                                      list(step = 67, finished = FALSE),
                                      NULL, NULL), 100),
    list(n = 179,
         tokens = list(bar_overall = "#++  ", p_running = "12% 67%"),
         result = FALSE))
  expect_equal(
    pmcmc_parallel_progress_data(list(list(step = 100, finished = TRUE),
                                      list(step = 100, finished = TRUE),
                                      list(step = 95, finished = FALSE)),
                                 100),
    list(n = 295,
         tokens = list(bar_overall = "##+", p_running = "95%"),
         result = FALSE))
  expect_equal(
    pmcmc_parallel_progress_data(list(list(step = 100, finished = TRUE),
                                      list(step = 100, finished = TRUE),
                                      list(step = 100, finished = TRUE)),
                                 100),
    list(n = 300,
         tokens = list(bar_overall = "###", p_running = ""),
         result = TRUE))
})


test_that("progress bar is a noop when progress = FALSE", {
  skip("rewrite")
  control <- pmcmc_control(100, n_workers = 2, n_chains = 3, progress = FALSE)
  p <- pmcmc_parallel_progress(control, force = TRUE)
  expect_silent(p(list(NULL, NULL, NULL)))
  expect_false(p(list(NULL, NULL, NULL)))
  expect_false(p(list(list(step = 100, finished = TRUE), NULL, NULL)))
  expect_silent(p(rep(list(list(step = 100, finished = TRUE)), 3)))
  expect_true(p(rep(list(list(step = 100, finished = TRUE)), 3)))
})


test_that("progress bar creates progress_bar when progress = TRUE", {
  skip("rewrite")
  control <- pmcmc_control(100, n_workers = 2, n_chains = 2, progress = TRUE)
  p <- pmcmc_parallel_progress(control, force = TRUE)

  Sys.sleep(0.2)
  expect_message(
    p(list(list(step = 50, finished = FALSE), NULL)),
    "\\[.\\] \\[\\+ \\] ETA .* \\| 00:00:[0-9]{2} so far \\(50%\\)")
  expect_message(
    p(list(list(step = 90, finished = FALSE),
           list(step = 10, finished = FALSE))),
    "\\[.\\] \\[\\+\\+\\] ETA .* \\| 00:00:[0-9]{2} so far \\(90% 10%\\)")
  expect_message(
    p(list(list(step = 100, finished = TRUE),
           list(step = 33, finished = FALSE))),
    "\\[.\\] \\[\\#\\+\\] ETA .* \\| 00:00:[0-9]{2} so far \\(33%\\)")
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
  control_parallel <- pmcmc_control(n_steps, n_chains = n_chains, n_workers = 2L)
  ans <- pmcmc(pars, filter, control = control_parallel)

  control_serial <- pmcmc_control(n_steps, n_chains = n_chains, progress = FALSE,
                                  use_parallel_seed = TRUE)
  cmp <- pmcmc(dat$pars, filter, control = control_serial)

  expect_equal(cmp$pars, ans$pars)
  expect_equal(cmp, ans)
})
