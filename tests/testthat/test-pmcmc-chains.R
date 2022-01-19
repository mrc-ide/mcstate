test_that("Split chain manually", {
  dat <- example_sir()
  n_particles <- 30
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  control <- pmcmc_control(10, n_chains = 4, use_parallel_seed = TRUE,
                           progress = FALSE)

  res1 <- pmcmc(dat$pars, p1, control = control)

  path <- tempfile()
  expect_equal(pmcmc_chains_prepare(path, dat$pars, p2, control),
               path)

  results <- vcapply(seq_len(control$n_chains), pmcmc_chains_run, path)
  expect_setequal(
    dir(path),
    c("control.rds", "inputs.rds", sprintf("results_%d.rds", 1:4)))

  res2 <- pmcmc_chains_collect(path)

  expect_equal(res1, res2)

  pmcmc_chains_cleanup(path)
  expect_false(file.exists(path))
})


test_that("can run split chains with nested model", {
  dat <- example_sir_shared()
  p1 <- particle_filter$new(dat$data, dat$model, 10, dat$compare,
                            dat$index, seed = 1L)
  p2 <- particle_filter$new(dat$data, dat$model, 10, dat$compare,
                            dat$index, seed = 1L)

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = matrix(0.00026),
    proposal_varied = matrix(0.00057))

  control <- pmcmc_control(10, n_chains = 2, use_parallel_seed = TRUE)
  res1 <- pmcmc(pars, p1, control = control)

  path <- tempfile()
  expect_equal(pmcmc_chains_prepare(path, dat$pars, p2, control),
               path)
  results <- vcapply(seq_len(control$n_chains), pmcmc_chains_run, path)
  res2 <- pmcmc_chains_collect(path)

  expect_equal(res1, res2)
})


test_that("split chain running requires one worker", {
  dat <- example_sir()
  n_particles <- 30
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  control <- pmcmc_control(10, n_chains = 4, n_workers = 2,
                           n_threads_total = 2,
                           use_parallel_seed = TRUE)
  expect_error(
    pmcmc_chains_prepare(tempfile(), dat$pars, p, control),
    "'n_workers' must be 1")
})


test_that("split chain running requires parallel seed setting", {
  dat <- example_sir()
  n_particles <- 30
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  control <- pmcmc_control(10, n_chains = 4, n_workers = 1,
                           use_parallel_seed = FALSE)
  expect_error(
    pmcmc_chains_prepare(tempfile(), dat$pars, p, control),
    "'use_parallel_seed' must be TRUE")
})


test_that("Error if results missing", {
  dat <- example_sir()
  n_particles <- 5
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  control <- pmcmc_control(3, n_chains = 4, use_parallel_seed = TRUE)
  path <- pmcmc_chains_prepare(tempfile(), dat$pars, p, control)
  expect_error(
    pmcmc_chains_collect(path),
    "Results missing for chains 1, 2, 3, 4")

  pmcmc_chains_run(2, path)
  pmcmc_chains_run(3, path)
  expect_error(
    pmcmc_chains_collect(path),
    "Results missing for chains 1, 4")
})


test_that("split chain running validates the chain id", {
  dat <- example_sir()
  n_particles <- 30
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  control <- pmcmc_control(10, n_chains = 4, n_workers = 1,
                           use_parallel_seed = TRUE)

  path <- pmcmc_chains_prepare(tempfile(), dat$pars, p, control)
  expect_error(
    pmcmc_chains_run(0, path),
    "'chain_id' must be at least 1",
    fixed = TRUE)
  expect_error(
    pmcmc_chains_run(5, path),
    "'chain_id' must be an integer in 1..4",
    fixed = TRUE)
})
