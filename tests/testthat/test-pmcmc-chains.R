test_that("Split chain manually", {
  dat <- example_sir()
  n_particles <- 30
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  control <- pmcmc_control(10, n_chains = 4, use_parallel_seed = TRUE)

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
})


test_that("error", {
  dat <- example_sir()
  n_particles <- 42
  n_steps <- 30
  n_chains <- 3

  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  control <- pmcmc_control(n_steps, n_chains = n_chains, n_workers = 2L,
                           n_steps_each = 20L)
  ## control$n_steps_each <- control$n_steps
  ## control <- pmcmc_control(n_steps, n_chains = n_chains)

  path <- tempfile()
  pmcmc_chains_prepare(path, dat$pars, p0, control, NULL)
  pmcmc_chains_run(1, path)
}




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
  skip("rework")
  dat <- example_sir()
  n_particles <- 30
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  control <- pmcmc_control(10, n_chains = 4, n_workers = 2,
                           use_parallel_seed = TRUE)
  expect_error(
    pmcmc_chains_prepare(dat$pars, p, NULL, control),
    "'n_workers' must be 1")
})


test_that("split chain running requires parallel seed setting", {
  skip("rework")
  dat <- example_sir()
  n_particles <- 30
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  control <- pmcmc_control(10, n_chains = 4, n_workers = 1,
                           use_parallel_seed = FALSE)
  expect_error(
    pmcmc_chains_prepare(dat$pars, p, NULL, control),
    "'use_parallel_seed' must be TRUE")
})


test_that("split chain running validates the chain id", {
  skip("rework")
  dat <- example_sir()
  n_particles <- 30
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  control <- pmcmc_control(10, n_chains = 4, n_workers = 1,
                           use_parallel_seed = TRUE)

  inputs <- pmcmc_chains_prepare(dat$pars, p, NULL, control)
  expect_error(
    pmcmc_chains_run(0, inputs),
    "'chain_id' must be at least 1",
    fixed = TRUE)
  expect_error(
    pmcmc_chains_run(5, inputs),
    "'chain_id' must be an integer in 1..4",
    fixed = TRUE)
})


test_that("Split chain and write to file", {
  skip("rework")
  dat <- example_sir()
  n_particles <- 30
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  control <- pmcmc_control(10, n_chains = 4, use_parallel_seed = TRUE)

  res1 <- pmcmc(dat$pars, p1, control = control)

  inputs <- pmcmc_chains_prepare(dat$pars, p2, NULL, control)
  ## typically the user will create this before but we won't here to
  ## show that this is robust to that
  path <- tempfile()

  samples_path <- vcapply(seq_len(control$n_chains), pmcmc_chains_run,
                          inputs, path)
  expect_true(file.exists(path))
  expect_true(file.info(path)$isdir)
  expect_true(all(file.exists(samples_path)))
  expect_setequal(dir(path), basename(samples_path))
  expect_equal(basename(samples_path),
               sprintf("samples_%d.rds", 1:4))
  samples <- lapply(samples_path, readRDS)
  res2 <- pmcmc_combine(samples = samples)

  expect_equal(res1, res2)
})
