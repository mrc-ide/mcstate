test_that("pmcmc_check_initial_nested - silent null initial", {
  dat <- example_uniform_shared()
  expect_equal(
    pmcmc_check_initial_nested(NULL, dat$pars, 2),
    array(dat$pars$initial(), c(4, 3, 2),
          dimnames = list(letters[1:4], c("p1", "p2", "p3"))))
})


test_that("pmcmc_check_initial_nested - silent matrix initial", {
  dat <- example_uniform_shared()
  expect_equal(
    pmcmc_check_initial_nested(matrix(1, 4, 3), dat$pars, 2),
    array(1, c(4, 3, 2),
          dimnames = list(letters[1:4], c("p1", "p2", "p3"), NULL)))
})


test_that("pmcmc_check_initial_nested - silent array initial", {
  dat <- example_uniform_shared()
  expect_equal(
    pmcmc_check_initial_nested(array(1, c(4, 3, 2)), dat$pars, 2),
    array(1, c(4, 3, 2),
          dimnames = list(letters[1:4], c("p1", "p2", "p3"), NULL)))
})


test_that("pmcmc_check_initial_nested - error array initial", {
  dat <- example_uniform_shared()
  expect_error(
    pmcmc_check_initial_nested(array(dim = c(4, 3, 1)), dat$pars, 2),
    "Expected 'initial' to be an array with dimensions 4 x 3 x 2")
  expect_error(
    pmcmc_check_initial_nested(array(dim = c(5, 3, 2)), dat$pars, 2),
    "Expected 'initial' to be an array with dimensions 4 x 3 x 2")
  expect_error(
    pmcmc_check_initial_nested(array(dim = c(4, 2, 2)), dat$pars, 2),
    "Expected 'initial' to be an array with dimensions 4 x 3 x 2")

  expect_error(
    pmcmc_check_initial_nested(
      array(dim = c(4, 3, 2),
            dimnames = list(letters[10:13], paste0("p", 1:3), NULL)),
      dat$pars, 2),
    "Expected names of dimension 1 of 'initial' to match parameters")
  expect_error(
    pmcmc_check_initial_nested(
      array(dim = c(4, 3, 2),
            dimnames = list(letters[1:4], letters[1:3], NULL)), dat$pars, 2),
    "Expected names of dimension 2 of 'initial' to match populations")
  expect_error(
    pmcmc_check_initial_nested(
      array(dim = c(4, 3, 2),
            dimnames = list(letters[1:4], NULL, letters[1:2])), dat$pars, 2),
    "Expected names of dimension 3 of 'initial' to be empty")

  expect_error(
    pmcmc_check_initial_nested(array(2, dim = c(4, 3, 2)), dat$pars, 2),
    "Starting point does not have finite prior probability (chain 1, 2)",
    fixed = TRUE)
})


test_that("pmcmc_check_initial_nested - error matrix initial", {
  dat <- example_uniform_shared()
  expect_error(
    pmcmc_check_initial_nested(matrix(0, 5, 3), dat$pars, 2),
    "Expected 'initial' to be a matrix with dimensions 4 x 3")
  expect_error(
    pmcmc_check_initial_nested(matrix(0, 4, 2), dat$pars, 2),
    "Expected 'initial' to be a matrix with dimensions 4 x 3")

  expect_error(
    pmcmc_check_initial_nested(
      matrix(0, 4, 3, dimnames = list(letters[2:5], letters[1:3])),
      dat$pars, 2),
    "Expected names of dimension 1 of 'initial' to match parameters")
  expect_error(
    pmcmc_check_initial_nested(
      matrix(0, 4, 3, dimnames = list(letters[1:4], letters[1:3])),
      dat$pars, 2),
    "Expected names of dimension 2 of 'initial' to match populations")

  expect_error(
    pmcmc_check_initial_nested(matrix(2, 4, 3), dat$pars, 2),
    "Starting point does not have finite prior probability (chain 1, 2)",
    fixed = TRUE)
})


test_that("pmcmc nested Uniform on unit square - fixed only", {
  testthat::skip_if_not_installed("coda")
  dat <- example_uniform_shared(varied = FALSE)
  control <- pmcmc_control(200, save_state = FALSE, save_trajectories = FALSE)

  set.seed(1)
  testthat::try_again(5, {
    res <- pmcmc(dat$pars, dat$filter, control = control)
    expect_s3_class(res, "mcstate_pmcmc")
    expect_true(all(acceptance_rate(res$pars[, , "p1"]) == 1))
    expect_true(all(acceptance_rate(res$pars[, , "p2"]) == 1))
    expect_true(all(acceptance_rate(res$pars[, , "p3"]) == 1))
    expect_true(abs(mean(res$pars[, "a", ]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars[, "b", ]) - 0.5) < 0.05)
  })
})


test_that("pmcmc nested Uniform on unit square - varied only", {
  testthat::skip_if_not_installed("coda")
  dat <- example_uniform_shared(fixed = FALSE)
  control <- pmcmc_control(200, save_state = FALSE, save_trajectories = FALSE)

  set.seed(1)
  testthat::try_again(5, {
    res <- pmcmc(dat$pars, dat$filter, control = control)
    expect_s3_class(res, "mcstate_pmcmc")
    expect_true(all(acceptance_rate(res$pars[, , "p1"]) == 1))
    expect_true(all(acceptance_rate(res$pars[, , "p2"]) == 1))
    expect_true(all(acceptance_rate(res$pars[, , "p3"]) == 1))
    expect_true(abs(mean(res$pars[, "c", ]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars[, "d", ]) - 0.5) < 0.05)
  })
})


test_that("pmcmc nested Uniform on unit square", {
  testthat::skip_if_not_installed("coda")
  dat <- example_uniform_shared()
  control <- pmcmc_control(201, save_state = FALSE, save_trajectories = FALSE)

  set.seed(1)
  testthat::try_again(5, {
    res <- pmcmc(dat$pars, dat$filter, control = control)
    expect_s3_class(res, "mcstate_pmcmc")
    expect_true(all(acceptance_rate(res$pars[, , "p1"]) == 0.5))
    expect_true(all(acceptance_rate(res$pars[, , "p2"]) == 0.5))
    expect_true(all(acceptance_rate(res$pars[, , "p3"]) == 0.5))
    expect_true(abs(mean(res$pars[, "a", ]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars[, "b", ]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars[, "c", ]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars[, "d", ]) - 0.5) < 0.05)
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
    expect_gt(ks_test(res$pars[i, "a", 1], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars[i, "a", 2], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars[i, "a", 3], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars[i, "b", 1], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars[i, "b", 2], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars[i, "b", 3], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars[i, "c", 1], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars[i, "c", 2], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars[i, "c", 3], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars[i, "d", 1], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars[i, "d", 2], "pnorm")$p.value, 0.05)
    expect_gt(ks_test(res$pars[i, "d", 3], "pnorm")$p.value, 0.05)
    expect_lt(abs(cov(res$pars[, , 1])[1, 2]), 0.1)
    expect_lt(abs(cov(res$pars[, , 2])[1, 2]), 0.1)
    expect_lt(abs(cov(res$pars[, , 3])[1, 2]), 0.1)
  })
})


test_that("nestedness must agree", {
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, 100, dat$compare,
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

  control <- pmcmc_control(5)
  expect_error(
    pmcmc(pars, p, control = control),
    "'pars' and 'filter' disagree on nestedness")
})


## We've ended up ~8% slower here for some of the abstraction.  It
## will be worth chasing these up later, but once the model moves more
## into the main code we'll be faster.  Using the compiled particle
## filter would also be much faster I think (though we currently lose
## the early exit).
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

  expect_equal(res1$pars[1, , ], pars$initial())
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
  expect_equal(res3$chain, rep(1:3, each = 50))
  expect_equal(res3$iteration, rep(1:50, 3))
  expect_equal(dim(res3$trajectories$state), c(3, 2, 150, 101))

  expect_equal(res1$pars, res3$pars[1:50, , ])
  expect_equal(res1$state, res3$state[, , 1:50])
  expect_equal(res1$restart$state, res3$restart$state[, , 1:50, ])
  expect_equal(res1$trajectories$state, res3$trajectories$state[, , 1:50, ])
})


test_that("return names on nested pmcmc output", {
  dat <- example_sir_shared()

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = matrix(0.00026),
    proposal_varied = matrix(0.00057))

  n_particles <- 10
  index2 <- function(info) list(run = 5L, state = c(a = 1, b = 2, c = 3))
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = index2, seed = 1)

  control <- pmcmc_control(5, save_state = TRUE, save_trajectories = TRUE)

  set.seed(1)
  results1 <- pmcmc(pars, p1, control = control)
  set.seed(1)
  results2 <- pmcmc(pars, p2, control = control)

  expect_null(rownames(results1$trajectories$state))
  expect_equal(rownames(results2$trajectories$state), c("a", "b", "c"))
  expect_equal(unname(results1$trajectories$state),
               unname(results2$trajectories$state))
})


test_that("run nested pmcmc with the particle filter and retain history", {
  dat <- example_sir_shared()

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = matrix(0.00026),
    proposal_varied = matrix(0.00057))


  n_particles <- 100
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)

  control1 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE)
  control2 <- pmcmc_control(30, save_trajectories = FALSE, save_state = FALSE)

  set.seed(1)
  results1 <- pmcmc(pars, p1, control = control1)
  set.seed(1)
  results2 <- pmcmc(pars, p2, control = control2)

  expect_setequal(
    names(results1),
    c("nested", "chain", "iteration", "pars_index",
      "pars", "probabilities", "state", "trajectories", "restart", "predict"))

  expect_null(results1$chain)
  expect_equal(results1$iteration, 1:30)

  ## Including or not the history does not change the mcmc trajectory:
  expect_identical(names(results1), names(results2))
  expect_equal(results1$pars, results2$pars)
  expect_equal(results1$probabilities, results2$probabilities)

  ## Parameters and probabilities have the expected shape
  expect_equal(dim(results1$pars), c(30, 2, 2))
  expect_equal(dimnames(results1$pars),
               list(NULL, c("beta", "gamma"), c("a", "b")))

  expect_equal(dim(results1$probabilities), c(30, 3, 2))
  expect_equal(
    dimnames(results1$probabilities),
    list(NULL, c("log_prior", "log_likelihood", "log_posterior"), c("a", "b")))

  ## History, if returned, has the correct shape
  expect_equal(dim(results1$state), c(5, 2, 30)) # state, pop, mcmc

  ## Trajectories, if returned, have the same shape
  expect_s3_class(results1$trajectories, "mcstate_trajectories")
  expect_equal(dim(results1$trajectories$state), c(3, 2, 30, 101))
  expect_equal(results1$trajectories$step, seq(0, 400, by = 4))
  expect_equal(results1$trajectories$rate, 4)

  ## Additional information required to predict
  expect_setequal(
    names(results1$predict),
    c("is_continuous", "transform", "index", "rate", "step", "time", "filter"))
})


test_that("nested_step_ratio works", {
  dat <- example_sir_shared()
  p <- particle_filter$new(dat$data, dat$model, 10, dat$compare,
                           dat$index)
  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = matrix(0.00026),
    proposal_varied = matrix(0.00057))

  ## Here, we never update beta, which is varied
  control <- pmcmc_control(30, nested_step_ratio = 30)
  res1 <- pmcmc(pars, p, control = control)
  expect_equal(as.numeric(res1$pars[, "beta", ]), rep(c(0.2, 0.3), each = 30))
  expect_equal(res1$pars[, "gamma", "a"], res1$pars[, "gamma", "b"])
  expect_false(all(res1$pars[, "gamma", "a"] == 0.1))

  ## Here, we never update gamma, which is fixed
  control <- pmcmc_control(30, nested_step_ratio = 1 / 30)
  res2 <- pmcmc(pars, p, control = control)
  expect_equal(res2$pars[, "gamma", ],
               matrix(0.1, 30, 2, dimnames = list(NULL, c("a", "b"))))
  expect_false(all(res2$pars[, "beta", ] == rep(c(0.2, 0.3), each = 30)))
})


test_that("Can do early exit with nested model", {
  dat <- example_sir_shared()
  n_particles <- 5
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            dat$index, seed = 1L)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            dat$index, seed = 1L)

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = matrix(0.00026),
    proposal_varied = matrix(0.00057))

  control1 <- pmcmc_control(30, filter_early_exit = FALSE)
  control2 <- pmcmc_control(30, filter_early_exit = TRUE)

  set.seed(1)
  res1 <- pmcmc(pars, p1, control = control1)
  set.seed(1)
  res2 <- pmcmc(pars, p2, control = control2)

  ## See similar test in test-pmcmc.R; we need to come up with a
  ## better test of this really.
  expect_false(identical(p2$inputs()$seed, p1$inputs()$seed))
})


test_that("Can run both fixed and varied at once", {
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

  control <- pmcmc_control(30, nested_step_ratio = 30,
                           nested_update_both = TRUE)
  res <- pmcmc(pars, p, control = control)

  ## Test that either all parameters change for everything or no
  ## parameters change everywhere:
  i <- which(res$pars[1:29, 1, 1] != res$pars[2:30, 1, 1])
  expect_true(all(res$pars[i, , ] != res$pars[i + 1, , ]))
  j <- setdiff(1:29, i)
  expect_true(all(res$pars[j, , ] == res$pars[j + 1, , ]))
})
