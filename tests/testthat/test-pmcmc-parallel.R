context("pmcmc (parallel)")

test_that("basic parallel operation", {
  dat <- example_sir()
  n_particles <- 42
  n_steps <- 30
  n_chains <- 3

  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  control <- pmcmc_control(n_steps, n_chains = n_chains, n_workers = 2L,
                           n_steps_each = 20L)
  ans <- pmcmc(dat$pars, p0, control = control)

  ## Run two chains manually with a given pair of seeds:
  s <- make_seeds(n_chains, 1L, dat$model)
  f <- function(idx) {
    set.seed(s[[idx]]$r)
    p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                             index = dat$index, seed = s[[idx]]$dust)
    pmcmc(dat$pars, p, control = pmcmc_control(n_steps))
  }

  samples <- lapply(seq_along(s), f)
  cmp <- pmcmc_combine(samples = samples)

  expect_equal(cmp$pars, ans$pars)
  expect_equal(cmp, ans)
})


test_that("Share out cores", {
  skip_on_cran()
  dat <- example_sir()
  n_particles <- 42
  n_steps <- 30
  n_chains <- 3

  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  control1 <- pmcmc_control(n_steps, n_chains = n_chains, n_workers = 2L,
                            n_steps_each = 15, n_threads_total = 4)
  ans <- pmcmc(dat$pars, p0, control = control1)

  control2 <- pmcmc_control(n_steps, n_chains = n_chains, n_workers = 1L,
                            n_steps_each = 15, use_parallel_seed = TRUE)
  cmp <- pmcmc(dat$pars, p0, control = control2)
  expect_equal(cmp$pars, ans$pars)
})


test_that("throw from callr operation", {
  skip_on_cran()

  dat <- example_sir()
  n_particles <- 42
  n_steps <- 30
  n_chains <- 3
  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  control <- pmcmc_control(n_steps, n_chains = 3, n_workers = 2,
                           n_steps_each = 5)

  initial <- pmcmc_check_initial(NULL, dat$pars, n_chains)
  seed <- make_seeds(n_chains, NULL, dat$model)

  control$n_workers <- 0
  obj <- pmcmc_orchestrator$new(dat$pars, initial, p0, control)
  path <- r6_private(obj)$path

  inputs <- readRDS(path$input)
  inputs$filter$n_threads <- "one"
  suppressWarnings(saveRDS(inputs, path$input))

  r <- pmcmc_remote$new(path$input, 2)
  r$wait_session_ready()
  r$init(1L)
  for (i in 1:20) {
    if (r$session$poll_process(1000) == "ready") {
      break
    } else {
      Sys.sleep(0.1)
    }
  }
  expect_error(r$read(), "'n_threads' must be an integer")
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


test_that("noop operations with a null thread pool", {
  skip_if_not_installed("mockery")
  mock_remote <- list2env(
    list(n_threads = 5,
         set_n_threads = mockery::mock()))

  p <- thread_pool$new(NULL, 2)
  expect_false(p$active)
  expect_equal(p$target, 1L)
  expect_null(p$add(10))
  expect_identical(p$free, 0L)
  expect_null(p$remove(mock_remote))
  expect_identical(p$free, 0L)
  mockery::expect_called(mock_remote$set_n_threads, 0L)
})


test_that("add and remove from a thread pool", {
  skip_if_not_installed("mockery")
  p <- thread_pool$new(20, 4)
  expect_true(p$active)
  expect_equal(p$n_workers, 4)
  expect_equal(p$n_threads, 20)
  expect_equal(p$free, 0)
  expect_equal(p$target, 5)

  mock_remote <- list2env(
    list(n_threads = 5,
         set_n_threads = mockery::mock()))

  p$add(mock_remote)
  expect_equal(p$free, 5)
  expect_equal(p$target, 7) # > 7 * 3 == 21

  p$remove(mock_remote)
  expect_equal(p$free, 3)
  p$remove(mock_remote)
  expect_equal(p$free, 1)
  p$remove(mock_remote)
  expect_equal(p$free, 0)
  p$remove(mock_remote)
  expect_equal(p$free, 0)

  mockery::expect_called(mock_remote$set_n_threads, 3L)
  expect_equal(
    mockery::mock_args(mock_remote$set_n_threads),
    list(list(7), list(7), list(6)))
})


test_that("construct parallel filter data", {
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
  control <- pmcmc_control(100, n_workers = 2, n_chains = 3, progress = FALSE)
  p <- pmcmc_parallel_progress(control, force = TRUE)
  expect_silent(p(list(NULL, NULL, NULL)))
  expect_false(p(list(NULL, NULL, NULL)))
  expect_false(p(list(list(step = 100, finished = TRUE), NULL, NULL)))
  expect_silent(p(rep(list(list(step = 100, finished = TRUE)), 3)))
  expect_true(p(rep(list(list(step = 100, finished = TRUE)), 3)))
})


test_that("progress bar creates progress_bar when progress = TRUE", {
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
  n_chains <- 3

  proposal_fixed <- matrix(0.00026)
  proposal_varied <- matrix(0.00057)

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = proposal_fixed, proposal_varied = proposal_varied)

  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  control <- pmcmc_control(n_steps, n_chains = n_chains, n_workers = 2L,
                           n_steps_each = 20L)
  ans <- pmcmc(pars, p0, control = control)

  ## Run two chains manually with a given pair of seeds:
  s <- make_seeds(n_chains, 1L, dat$model)
  f <- function(idx) {
    set.seed(s[[idx]]$r)
    p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                             index = dat$index, seed = s[[idx]]$dust)
    pmcmc(dat$pars, p, control = pmcmc_control(n_steps))
  }

  samples <- lapply(seq_along(s), f)
  cmp <- pmcmc_combine(samples = samples)

  expect_equal(cmp$pars, ans$pars)
  expect_equal(cmp$state, ans$state)
})


## Recent versions of callr/processx/covr have conspired to not reveal
## this bit of logic through an integration test
test_that("throw from callr operation", {
  skip_on_cran()

  dat <- example_sir()
  n_particles <- 42
  n_steps <- 30
  n_chains <- 3
  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  control <- pmcmc_control(n_steps, n_chains = 3, n_workers = 2,
                           n_steps_each = 5)

  initial <- pmcmc_check_initial(NULL, dat$pars, n_chains)
  seed <- make_seeds(n_chains, NULL, dat$model)

  control$n_workers <- 0
  obj <- pmcmc_orchestrator$new(dat$pars, initial, p0, control)
  path <- r6_private(obj)$path

  r <- pmcmc_remote$new(path$input, 1)
  r$wait_session_ready()
  r$init(1L)
  for (i in 1:20) {
    if (r$session$poll_process(1000) == "ready") {
      break
    } else {
      Sys.sleep(0.1)
    }
  }
  r$read()
  expect_equal(r$n_threads, 1)

  prev <- r$set_n_threads(2)
  expect_equal(prev, 1)
  expect_equal(r$n_threads, 2)

  prev <- r$set_n_threads(1)
  expect_equal(prev, 2)
  expect_equal(r$n_threads, 1)

  r$session$kill(1)
})
