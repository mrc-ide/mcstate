test_that("pmcmc_check_initial_nested - silent null initial", {
  dat <- example_uniform_shared()
  expect_equal(
    pmcmc_check_initial_nested(NULL, dat$pars, 2),
    array(dat$pars$initial(), c(3, 4, 2),
          dimnames = list(c("p1", "p2", "p3"), letters[1:4]))
  )
})

test_that("pmcmc_check_initial_nested - silent matrix initial", {
  dat <- example_uniform_shared()
  expect_equal(
    pmcmc_check_initial_nested(matrix(1, 3, 4), dat$pars, 2),
    array(1, c(3, 4, 2), dimnames = list(paste0("p", 1:3), letters[1:4], NULL))
  )
})

test_that("pmcmc_check_initial_nested - silent array initial", {
  dat <- example_uniform_shared()
  expect_equal(
    pmcmc_check_initial_nested(array(1, c(3, 4, 2)), dat$pars, 2),
    array(1, c(3, 4, 2), dimnames = list(paste0("p", 1:3), letters[1:4], NULL))
  )
})

test_that("pmcmc_check_initial_nested - error array initial", {
  dat <- example_uniform_shared()
  expect_error(
    pmcmc_check_initial_nested(array(dim = c(3, 4, 1)), dat$pars, 2),
    "Expected an array with 2 layers"
  )
  expect_error(
    pmcmc_check_initial_nested(array(dim = c(3, 5, 2)), dat$pars, 2),
    "Expected an array with 4 columns"
  )
  expect_error(
    pmcmc_check_initial_nested(array(dim = c(2, 4, 2)), dat$pars, 2),
    "Expected an array with 3 rows"
  )
  expect_error(
    pmcmc_check_initial_nested(
      array(dim = c(3, 4, 2),
            dimnames = list(letters[1:3], letters[1:4], NULL)), dat$pars, 2),
    "has rownames"
  )
  expect_error(
    pmcmc_check_initial_nested(
      array(dim = c(3, 4, 2),
            dimnames = list(paste0("p", 1:3), letters[10:13], NULL)),
      dat$pars, 2),
    "has colnames"
  )
  expect_error(
    pmcmc_check_initial_nested(array(2, dim = c(3, 4, 2)), dat$pars, 2),
    "finite prior"
  )
})

test_that("pmcmc_check_initial_nested - error matrix initial", {
  dat <- example_uniform_shared()
  expect_error(
    pmcmc_check_initial_nested(matrix(nrow = 3, ncol = 5), dat$pars, 2),
    "Expected a matrix with 4 columns"
  )
  expect_error(
    pmcmc_check_initial_nested(matrix(nrow = 2, ncol = 4), dat$pars, 2),
    "Expected a matrix with 3 rows"
  )
  expect_error(
    pmcmc_check_initial_nested(
      matrix(nrow = 3, ncol = 4,
             dimnames = list(letters[1:3], letters[1:4])), dat$pars, 2),
    "has rownames"
  )
  expect_error(
    pmcmc_check_initial_nested(
      matrix(nrow = 3, ncol = 4,
             dimnames = list(paste0("p", 1:3), letters[10:13])), dat$pars, 2),
    "has colnames"
  )
  expect_error(pmcmc_check_initial_nested(matrix(2, 3, 4), dat$pars, 2),
               "finite prior")
})


test_that("pmcmc nested Uniform on unit square - fixed only", {
  dat <- example_uniform_shared(varied = FALSE)
  control <- pmcmc_control(1000, save_state = FALSE, save_trajectories = FALSE)

  set.seed(1)
  testthat::try_again(5, {
    res <- pmcmc(dat$pars, dat$filter, control = control)
    expect_s3_class(res, "mcstate_pmcmc")
    expect_true(all(acceptance_rate(t(res$pars[, "p1", ])) == 0.5))
    expect_true(all(acceptance_rate(t(res$pars[, "p2", ])) == 0.5))
    expect_true(all(acceptance_rate(t(res$pars[, "p3", ])) == 0.5))
    expect_true(abs(mean(res$pars["a", , ]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars["b", , ]) - 0.5) < 0.05)
  })
})

test_that("pmcmc nested Uniform on unit square - varied only", {
  dat <- example_uniform_shared(fixed = FALSE)
  control <- pmcmc_control(1000, save_state = FALSE, save_trajectories = FALSE)

  set.seed(1)
  testthat::try_again(5, {
    res <- pmcmc(dat$pars, dat$filter, control = control)
    expect_s3_class(res, "mcstate_pmcmc")
    expect_true(all(acceptance_rate(t(res$pars[, "p1", ])) == 0.5))
    expect_true(all(acceptance_rate(t(res$pars[, "p2", ])) == 0.5))
    expect_true(all(acceptance_rate(t(res$pars[, "p3", ])) == 0.5))
    expect_true(abs(mean(res$pars["c", , ]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars["d", , ]) - 0.5) < 0.05)
  })
})


test_that("pmcmc nested Uniform on unit square", {
  dat <- example_uniform_shared()
  control <- pmcmc_control(1000, save_state = FALSE, save_trajectories = FALSE)

  set.seed(1)
  testthat::try_again(5, {
    res <- pmcmc(dat$pars, dat$filter, control = control)
    expect_s3_class(res, "mcstate_pmcmc")
    expect_true(all(acceptance_rate(t(res$pars[, "p1", ])) == 0.5))
    expect_true(all(acceptance_rate(t(res$pars[, "p2", ])) == 0.5))
    expect_true(all(acceptance_rate(t(res$pars[, "p3", ])) == 0.5))
    expect_true(abs(mean(res$pars["a", , ]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars["b", , ]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars["c", , ]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars["d", , ]) - 0.5) < 0.05)
  })
})

test_that("pmcmc nested multivariate gaussian", {
  testthat::skip_if_not_installed("mvtnorm")
  dat <- example_mvnorm_shared()
  control <- pmcmc_control(500, save_state = FALSE, save_trajectories = FALSE)

  set.seed(1)
  testthat::try_again(10, {
    res <- pmcmc(dat$pars, dat$filter, control = control)
    i <- seq(1, 500, by = 20)
    expect_s3_class(res, "mcstate_pmcmc")
    ks_test <- function(x, y) suppressWarnings(ks.test(x, y))
    expect_gt(ks_test(res$pars["a", 1, i], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars["a", 2, i], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars["a", 3, i], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars["b", 1, i], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars["b", 2, i], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars["b", 3, i], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars["c", 1, i], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars["c", 2, i], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars["c", 3, i], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars["d", 1, i], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars["d", 2, i], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars["d", 3, i], "pnorm")$p.value, 0.05)
    expect_lt(abs(cov(res$pars[, 1, ])[1, 2]), 0.1)
    expect_lt(abs(cov(res$pars[, 2, ])[1, 2]), 0.1)
    expect_lt(abs(cov(res$pars[, 3, ])[1, 2]), 0.1)
  })
})

test_that("pmcmc nested sir - 1 chain", {
  dat <- example_sir_shared()
  p <- particle_filter$new(dat$data, dat$model, 100, dat$compare,
                           dat$index)
  control <- pmcmc_control(30, save_state = TRUE, save_trajectories = TRUE,
                           save_restart = TRUE, rerun_every = 10,
                           n_threads_total = 2, use_parallel_seed = TRUE,
                           n_chains = 1L)
  proposal_fixed <- matrix(0.00026)
  proposal_varied <- matrix(0.00057)

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = proposal_fixed, proposal_varied = proposal_varied)

  set.seed(1)
  res1 <- pmcmc(pars, p, control = control)
  set.seed(1)
  res2 <- pmcmc(pars, p, control = control)
  expect_equal(res1, res2)
})


test_that("pmcmc nested sir - 2 chains", {
  dat <- example_sir_shared()
  p1 <- particle_filter$new(dat$data, dat$model, 100, dat$compare,
                            dat$index, seed = 1L)
  p2 <- particle_filter$new(dat$data, dat$model, 100, dat$compare,
                            dat$index, seed = 1L)
  proposal_fixed <- matrix(0.00026)
  proposal_varied <- matrix(0.00057)

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = proposal_fixed, proposal_varied = proposal_varied)

  control1 <- pmcmc_control(50, save_state = TRUE, n_chains = 1,
                            save_restart = c(10, 20, 30, 40),
                            save_trajectories = TRUE)
  control2 <- pmcmc_control(50, n_chains = 3, save_state = TRUE,
                            save_restart = c(10, 20, 30, 40),
                            save_trajectories = TRUE)

  set.seed(1)
  res1 <- pmcmc(pars, p1, control = control1)
  expect_s3_class(res1, "mcstate_pmcmc")
  expect_null(res1$chain)

  set.seed(1)
  res3 <- pmcmc(pars, p2, control = control2)
  expect_s3_class(res3, "mcstate_pmcmc")
  expect_equal(res3$chain, rep(1:3, each = 51))
  expect_equal(dim(res3$trajectories$state), c(3, 153, 2, 101))

  expect_equal(res1$pars, res3$pars[, , 1:51])
  expect_equal(res1$state, res3$state[, , 1:51])
  expect_equal(res1$restart$state, res3$restart$state[, 1:51, , ])
  expect_equal(res1$trajectories$state, res3$trajectories$state[, 1:51, , ])
})

test_that("sample_trajectory_nested single state", {
  expect_equal(sample_trajectory_nested(array(1, c(1, 1, 1, 1))),
               array(1, c(1, 1, 1)))
})


test_that("return names on nested pmcmc output", {
  dat <- example_sir_shared()
  p1 <- particle_filter$new(dat$data, dat$model, 100, dat$compare,
                            dat$index, seed = 1L)

  proposal_fixed <- matrix(0.00026)
  proposal_varied <- matrix(0.00057)

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = proposal_fixed, proposal_varied = proposal_varied)

  control1 <- pmcmc_control(50, save_state = FALSE, n_chains = 1)

  n_particles <- 10
  index2 <- function(info) list(run = 4L, state = c(a = 1, b = 2, c = 3))
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = index2)

  control <- pmcmc_control(5, save_state = TRUE, save_trajectories = TRUE)

  set.seed(1)
  results1 <- pmcmc(pars, p1, control = control)
  set.seed(1)
  results2 <- pmcmc(pars, p2, control = control)

  expect_null(rownames(results1$trajectories$state))
  expect_equal(rownames(results2$trajectories$state), c("a", "b", "c"))
})



test_that("run nested pmcmc with the particle filter and retain history", {
  dat <- example_sir_shared()

  proposal_fixed <- matrix(0.00026)
  proposal_varied <- matrix(0.00057)

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = proposal_fixed, proposal_varied = proposal_varied)


  n_particles <- 100
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)

  control1 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE)
  control2 <- pmcmc_control(30, save_trajectories = FALSE, save_state = FALSE)

  set.seed(1)
  results1 <- pmcmc(pars, p1, control = control1)
  set.seed(1)
  results2 <- pmcmc(pars, p2, control = control2)

  expect_setequal(
    names(results1),
    c("chain", "iteration",
      "pars", "probabilities", "state", "trajectories", "restart", "predict"))

  expect_null(results1$chain)
  expect_equal(results1$iteration, 0:30)

  ## Including or not the history does not change the mcmc trajectory:
  expect_identical(names(results1), names(results2))
  expect_equal(results1$pars, results2$pars)
  expect_equal(results1$probabilities, results2$probabilities)

  ## Parameters and probabilities have the expected shape
  expect_equal(dim(results1$pars), c(2, 2, 31))
  expect_equal(rownames(results1$pars), c("beta", "gamma"))
  expect_equal(colnames(results1$pars), c("a", "b"))

  expect_equal(dim(results1$probabilities), c(3, 2, 31))
  expect_equal(rownames(results1$probabilities),
               c("log_prior", "log_likelihood", "log_posterior"))

  ## History, if returned, has the correct shape
  expect_equal(dim(results1$state), c(5, 2, 31)) # state, mcmc

  ## Trajectories, if returned, have the same shape
  expect_s3_class(results1$trajectories, "mcstate_trajectories")
  expect_equal(dim(results1$trajectories$state), c(3, 31, 2, 101))
  expect_equal(results1$trajectories$step, seq(0, 400, by = 4))
  expect_equal(results1$trajectories$rate, 4)

  ## Additional information required to predict
  expect_setequal(
    names(results1$predict),
    c("transform", "index", "rate", "step", "filter"))
})


test_that("nested_step_ratio works", {
  dat <- example_sir_shared()
  p <- particle_filter$new(dat$data, dat$model, 10, dat$compare,
                           dat$index)
  proposal_fixed <- matrix(0.00026)
  proposal_varied <- matrix(0.00057)

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = proposal_fixed, proposal_varied = proposal_varied)

  control <- pmcmc_control(30, nested_step_ratio = 30)
  res1 <- pmcmc(pars, p, control = control)
  expect_equal(as.numeric(res1$pars[1, , ]), rep(c(0.2, 0.3), 31))
  expect_false(identical(as.numeric(res1$pars[2, , ]), rep(0.1, 62)))

  control <- pmcmc_control(30, nested_step_ratio = 1 / 30)
  res2 <- pmcmc(pars, p, control = control)
  expect_equal(as.numeric(res2$pars[2, , ]), rep(0.1, 62))
  expect_false(identical(as.numeric(res2$pars[1, , ]), rep(c(0.2, 0.3), 31)))

  control <- pmcmc_control(30, nested_step_ratio = 1 / 2)
  expect_is(pmcmc(pars, p, control = control), "mcstate_pmcmc")

  control <- pmcmc_control(30, nested_step_ratio = 2)
  expect_is(pmcmc(pars, p, control = control), "mcstate_pmcmc")
})


test_that("can run split chains with nested model", {
  dat <- example_sir_shared()
  p1 <- particle_filter$new(dat$data, dat$model, 10, dat$compare,
                            dat$index, seed = 1L)
  p2 <- particle_filter$new(dat$data, dat$model, 10, dat$compare,
                            dat$index, seed = 1L)
  proposal_fixed <- matrix(0.00026)
  proposal_varied <- matrix(0.00057)

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = proposal_fixed, proposal_varied = proposal_varied)

  control <- pmcmc_control(10, n_chains = 2, use_parallel_seed = TRUE)
  res1 <- pmcmc(pars, p1, control = control)

  inputs <- pmcmc_chains_prepare(dat$pars, p2, NULL, control)
  samples <- lapply(seq_len(control$n_chains), pmcmc_chains_run, inputs)
  res2 <- pmcmc_combine(samples = samples)

  v <- c("chain", "iteration", "pars", "probabilities", "trajectories")
  expect_equal(res1[v], res2[v])
})
