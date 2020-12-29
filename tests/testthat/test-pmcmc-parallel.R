context("pmmc (parallel)")

test_that("basic parallel operation", {
  dat <- example_sir()
  n_particles <- 42
  n_steps <- 30
  n_chains <- 3

  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  obj <- orchestrator$new(dat$pars, p0, n_steps, n_chains = n_chains,
                          n_workers = 2, n_steps_each = 20)
  while (!obj$step()) {
  }
  ans <- obj$get_results()

  ## Run two chains manually with a given pair of seeds:
  s <- make_seeds(n_chains, 1L)
  f <- function(idx) {
    set.seed(s[[idx]]$r)
    p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                              index = dat$index, seed = s[[idx]]$dust)
    pmcmc(dat$pars, p, n_steps)
  }

  samples <- lapply(seq_along(s), f)
  cmp <- pmcmc_combine(samples = samples)

  expect_equal(cmp$pars, ans$pars)
  expect_equal(cmp, ans)
})


test_that("make seeds with integer returns length-32 raws", {
  s <- make_seeds(2, 1)
  expect_length(s, 2L)
  expect_setequal(names(s[[1]]), c("dust", "r"))
  expect_length(s[[1]]$dust, 32L)
  expect_is(s[[1]]$dust, "raw")
  expect_identical(make_seeds(2, 1), s)
  expect_identical(make_seeds(4, 1)[1:2], s)
})


test_that("make seeds with long raw retains size", {
  seed <- dust::dust_rng$new(NULL, n_generators = 3)$state()
  s <- make_seeds(2L, seed)
  expect_length(s, 2L)
  expect_setequal(names(s[[1]]), c("dust", "r"))
  expect_identical(s[[1]]$dust, seed)
  expect_length(s[[2]]$dust, length(seed))
})


test_that("Don't run with fewer chains than workers", {
  dat <- example_sir()
  n_particles <- 42
  n_steps <- 30
  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  expect_error(
    obj <- orchestrator$new(dat$pars, p0, n_steps, n_chains = 2,
                            n_workers = 5, n_steps_each = 20),
    "'n_chains' (2) is less than 'n_workers' (5)",
    fixed = TRUE)
})


test_that("Don't run with invalid n_threads", {
  dat <- example_sir()
  n_particles <- 42
  n_steps <- 30
  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  expect_error(
    obj <- orchestrator$new(dat$pars, p0, n_steps, n_chains = 5,
                            n_workers = 5, n_steps_each = 20, n_threads = 3),
    "'n_threads' (3) is less than 'n_workers' (5)",
    fixed = TRUE)
  expect_error(
    obj <- orchestrator$new(dat$pars, p0, n_steps, n_chains = 5,
                            n_workers = 5, n_steps_each = 20, n_threads = 8),
    "'n_threads' (8) is not a multiple of 'n_workers' (5)",
    fixed = TRUE)
})


test_that("Share out cores", {
  skip_on_cran()
  dat <- example_sir()
  n_particles <- 42
  n_steps <- 30
  n_chains <- 3

  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  obj <- orchestrator$new(dat$pars, p0, n_steps, n_chains = n_chains,
                          n_workers = 2, n_steps_each = 15, n_threads = 4)
  while (!obj$step()) {
  }

  ans <- obj$get_results()
})


test_that("noop operations with a null thread pool", {
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


test_that("throw from callr operation", {
  skip_on_cran()
  dat <- example_sir()
  n_particles <- 42
  n_steps <- 30
  n_chains <- 3
  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  inputs <- p0$inputs()
  control <- list(n_steps = 30, n_steps_each = 5, rerun_every = Inf,
                  save_state = FALSE, save_trajectories = FALSE)
  initial <- pmcmc_check_initial(NULL, dat$pars, n_chains)
  seed <- make_seeds(n_chains, NULL)
  inputs$n_threads <- "one"
  r <- remote$new(dat$pars, initial, inputs, control, seed)
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
