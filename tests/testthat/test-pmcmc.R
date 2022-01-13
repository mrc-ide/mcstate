context("pmcmc")


## TODO: Ed and Lilith we could use better tests throughout here. The
## sampler should be ok to run for ~10k iterations without taking too
## long to be annoying in tests.
test_that("mcmc works for uniform distribution on unit square", {
  dat <- example_uniform()
  control <- pmcmc_control(1000, save_state = FALSE, save_trajectories = FALSE)

  set.seed(1)
  testthat::try_again(5, {
    res <- pmcmc(dat$pars, dat$filter, control = control)
    expect_s3_class(res, "mcstate_pmcmc")
    expect_true(all(acceptance_rate(res$pars) == 1))
    expect_true(abs(mean(res$pars[, "a"]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars[, "b"]) - 0.5) < 0.05)
  })
})


test_that("mcmc works for multivariate gaussian", {
  testthat::skip_if_not_installed("mvtnorm")
  dat <- example_mvnorm()
  control <- pmcmc_control(1000, save_state = FALSE, save_trajectories = FALSE)

  set.seed(1)
  testthat::try_again(5, {
    res <- pmcmc(dat$pars, dat$filter, control = control)
    i <- seq(1, 1000, by = 20)
    expect_s3_class(res, "mcstate_pmcmc")
    expect_gt(ks.test(res$pars[i, "a"], "pnorm")$p.value, 0.05)
    expect_gt(ks.test(res$pars[i, "b"], "pnorm")$p.value, 0.05)
    expect_lt(abs(cov(res$pars)[1, 2]), 0.1)
  })
})


test_that("Reflect parameters: double sided", {
  expect_equal(reflect_proposal(0.4, 0, 1), 0.4)
  expect_equal(reflect_proposal(1.4, 0, 1), 0.6)
  expect_equal(reflect_proposal(2.4, 0, 1), 0.4)

  ## Additional cases from the original
  expect_equal(reflect_proposal(6, 1, 5), 4)
  expect_equal(reflect_proposal(0, 1, 5), 2)
  expect_equal(reflect_proposal(10, 1, 5), 2)

  expect_equal(reflect_proposal(0:2 + 0.4, rep(0, 3), rep(1, 3)),
               c(0.4, 0.6, 0.4))
  expect_equal(reflect_proposal(0:2 + 1.4, rep(1, 3), rep(2, 3)),
               c(1.4, 1.6, 1.4))
})


test_that("Reflect parameters: infinite", {
  p <- rnorm(10)
  expect_equal(reflect_proposal(p, rep(-Inf, 5), rep(Inf, 5)), p)
})


test_that("reflect parameters: lower", {
  expect_equal(reflect_proposal(-0.4, 0, Inf), 0.4)
  expect_equal(reflect_proposal(0.6, 1, Inf), 1.4)
  expect_equal(reflect_proposal(0.6, 0, Inf), 0.6)
})


test_that("reflect parameters: upper", {
  expect_equal(reflect_proposal(1.4, -Inf, 1), 0.6)
  expect_equal(reflect_proposal(0.4, -Inf, 0), -0.4)
  expect_equal(reflect_proposal(0.4, -Inf, 1), 0.4)
})


test_that("proposal uses provided covariance structure", {
  dat <- example_uniform(proposal_kernel = matrix(0.1, 2, 2))
  control <- pmcmc_control(100, save_state = FALSE, save_trajectories = FALSE)
  res <- pmcmc(dat$pars, dat$filter, control = control)

  expect_true(all(acceptance_rate(res$pars) == 1))
  expect_equal(res$pars[, 1], res$pars[, 2])
})


test_that("run pmcmc with the particle filter and retain history", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  pars <- pmcmc_parameters$new(
    list(pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                         prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal = proposal_kernel)

  dat <- example_sir()
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
    c("nested", "chain", "iteration", "pars_index",
      "pars", "probabilities", "state", "trajectories", "restart", "predict"))

  expect_null(results1$pars_index)
  expect_null(results1$chain)
  expect_equal(results1$iteration, 1:30)

  ## Including or not the history does not change the mcmc trajectory:
  expect_identical(names(results1), names(results2))
  expect_equal(results1$pars, results2$pars)
  expect_equal(results1$probabilities, results2$probabilities)

  ## We did mix
  expect_true(all(acceptance_rate(results1$pars) > 0))

  ## Parameters and probabilities have the expected shape
  expect_equal(dim(results1$pars), c(30, 2))
  expect_equal(colnames(results1$pars), c("beta", "gamma"))

  expect_equal(dim(results1$probabilities), c(30, 3))
  expect_equal(colnames(results1$probabilities),
               c("log_prior", "log_likelihood", "log_posterior"))

  ## History, if returned, has the correct shape
  expect_equal(dim(results1$state), c(5, 30)) # state, mcmc

  ## Trajectories, if returned, have the same shape
  expect_s3_class(results1$trajectories, "mcstate_trajectories")
  expect_equal(dim(results1$trajectories$state), c(3, 30, 101))
  expect_equal(
    results1$trajectories$state[, , dim(results1$trajectories$state)[3]],
    results1$state[1:3, ])
  expect_equal(results1$trajectories$predicted, rep(FALSE, 101))
  expect_equal(results1$trajectories$step, seq(0, 400, by = 4))
  expect_equal(results1$trajectories$rate, 4)

  ## Additional information required to predict
  expect_setequal(
    names(results1$predict),
    c("transform", "index", "rate", "step", "filter"))
  expect_identical(results1$predict$transform, as.list)
  expect_equal(results1$predict$index, 1:3)
  expect_equal(results1$predict$rate, 4)
  expect_equal(results1$predict$step, last(dat$data$step_end))
  expect_identical(results1$predict$filter, p1$inputs())
})


test_that("collecting state from model yields an RNG state", {
  dat <- example_sir()
  n_particles <- 30
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  p3 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 2L)

  set.seed(1)
  results1 <- pmcmc(dat$pars, p1,
                    control = pmcmc_control(5, save_state = FALSE))
  set.seed(1)
  results2 <- pmcmc(dat$pars, p2,
                    control = pmcmc_control(5, save_state = TRUE))
  set.seed(1)
  results3 <- pmcmc(dat$pars, p3,
                    control = pmcmc_control(5, save_state = TRUE))

  expect_identical(
    r6_private(p1)$last_model$rng_state(),
    r6_private(p2)$last_model$rng_state())
  expect_identical(
    results2$predict$filter$seed,
    r6_private(p1)$last_model$rng_state()[1:32])
  expect_false(
    identical(results2$predict$filter$seed, results3$predict$filter$seed))
})


test_that("running pmcmc with progress = TRUE prints messages", {
  dat <- example_uniform()
  control <- pmcmc_control(1000, save_state = FALSE, save_trajectories = FALSE,
                           progress = TRUE)
  expect_message(
    pmcmc(dat$pars, dat$filter, control = control),
    "Finished 1000 steps in ")
})


test_that("run multiple chains", {
  dat <- example_uniform()

  control1 <- pmcmc_control(100, save_state = FALSE, n_chains = 1)
  control2 <- pmcmc_control(100, save_state = FALSE, n_chains = 3)

  set.seed(1)
  res1 <- pmcmc(dat$pars, dat$filter, control = control1)
  expect_s3_class(res1, "mcstate_pmcmc")
  expect_null(res1$chain)

  set.seed(1)
  res3 <- pmcmc(dat$pars, dat$filter, control = control2)
  expect_s3_class(res3, "mcstate_pmcmc")
  expect_equal(res3$chain, rep(1:3, each = 100))

  expect_equal(res1$pars, res3$pars[1:100, ])
})


test_that("progress in multiple chains", {
  dat <- example_uniform()
  control <- pmcmc_control(100, save_state = FALSE, n_chains = 3,
                           progress = TRUE)
  expect_message(
    pmcmc(dat$pars, dat$filter, control = control),
    "Running chain 2 / 3")
})


test_that("return names on pmcmc output", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  pars <- pmcmc_parameters$new(
    list(pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                         prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal = proposal_kernel)

  dat <- example_sir()
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


test_that("can use default initial conditions", {
  dat <- example_sir()
  dn <- list(dat$pars$names(), NULL)
  expect_equal(pmcmc_check_initial(NULL, dat$pars, 1),
               matrix(c(0.2, 0.1), 2, 1, dimnames = dn))
  expect_equal(pmcmc_check_initial(NULL, dat$pars, 5),
               matrix(c(0.2, 0.1), 2, 5, dimnames = dn))
})


test_that("can use a vector initial conditions and expand it out", {
  dat <- example_uniform()
  dn <- list(dat$pars$names(), NULL)
  expect_equal(pmcmc_check_initial(c(0.1, 0.2), dat$pars, 1),
               matrix(c(0.1, 0.2), 2, 1, dimnames = dn))
  expect_equal(pmcmc_check_initial(c(0.1, 0.2), dat$pars, 5),
               matrix(c(0.1, 0.2), 2, 5, dimnames = dn))
  expect_equal(pmcmc_check_initial(c(a = 0.1, b = 0.2), dat$pars, 5),
               matrix(c(0.1, 0.2), 2, 5, dimnames = dn))
})


test_that("can validate a vector of initial conditions", {
  dat <- example_uniform()
  expect_error(pmcmc_check_initial(c(0.1, 0.2, 0.4), dat$pars, 1),
               "Expected 'initial' to be a vector with length 2")
  expect_error(pmcmc_check_initial(c(x = 0.1, y = 0.2), dat$pars, 1),
               "Expected names of 'initial' to match parameters ('a', 'b')",
               fixed = TRUE)
  expect_error(pmcmc_check_initial(c(-0.1, 0.2), dat$pars, 1),
               "Starting point does not have finite prior probability",
               fixed = TRUE)
})


test_that("can use a matrix initial conditions", {
  dat <- example_uniform()
  dn <- list(dat$pars$names(), NULL)
  expect_equal(pmcmc_check_initial(cbind(c(0.1, 0.2)), dat$pars, 1),
               matrix(c(0.1, 0.2), 2, 1, dimnames = dn))
  m <- matrix(runif(10), 2, 5, dimnames = dn)
  expect_equal(pmcmc_check_initial(unname(m), dat$pars, 5), m)
  expect_equal(pmcmc_check_initial(m, dat$pars, 5), m)
})


test_that("can validate a matrix initial conditions", {
  dat <- example_uniform()
  expect_error(
    pmcmc_check_initial(matrix(0.5, 3, 5), dat$pars, 5),
    "Expected 'initial' to be a matrix with dimensions 2 x 5")
  expect_error(
    pmcmc_check_initial(matrix(0.5, 2, 6), dat$pars, 5),
    "Expected 'initial' to be a matrix with dimensions 2 x 5")

  expect_error(
    pmcmc_check_initial(matrix(0.5, 2, 5, dimnames = list(c("x", "y"), NULL)),
                        dat$pars, 5),
    "Expected names of dimension 1 of 'initial' to match parameters ('a', 'b')",
    fixed = TRUE)

  m <- matrix(runif(10), 2, 5)
  i <- cbind(c(2, 1, 2), c(2, 4, 5))
  m[i] <- -m[i]
  expect_error(
    pmcmc_check_initial(m, dat$pars, 5),
    "Starting point does not have finite prior probability (chain 2, 4, 5)",
    fixed = TRUE)
})


test_that("can start a pmcmc from a matrix of starting points", {
  proposal_kernel <- diag(2) * 2
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")
  pars <- pmcmc_parameters$new(
    list(pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                         prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal = proposal_kernel * 50)

  p0 <- pars$initial()
  initial <- matrix(p0, 2, 3, dimnames = list(names(p0), NULL))
  initial[] <- initial + runif(6, 0, 0.001)

  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, 10, dat$compare,
                           index = dat$index)

  control <- pmcmc_control(2, n_chains = 3)
  res <- pmcmc(pars, p, control = control, initial = initial)

  expect_equal(nrow(unique(res$pars[res$iteration == 1, ])), 3)
})


test_that("can trigger rerunning particle filter", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  pars <- pmcmc_parameters$new(
    list(pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                         prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal = proposal_kernel * 50)

  dat <- example_sir()
  n_particles <- 100
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)

  control1 <- pmcmc_control(30, save_state = TRUE, save_trajectories = TRUE,
                            rerun_every = 2)
  control2 <- pmcmc_control(30, save_state = TRUE, save_trajectories = TRUE,
                            rerun_every = Inf)

  set.seed(1)
  res1 <- pmcmc(pars, p1, control = control1)
  expect_lte(max(rle(res1$probabilities[, "log_likelihood"])$lengths), 2)

  set.seed(1)
  res2 <- pmcmc(pars, p1, control = control2)
  expect_gt(max(rle(res2$probabilities[, "log_likelihood"])$lengths), 5)
})


test_that("rerunning the particle filter triggers the filter run method", {
  skip_if_not_installed("mockery")
  dat <- example_uniform()
  dat$filter$run <- mockery::mock(1, cycle = TRUE)
  dat$inputs <- function() NULL

  ## with the dummy version we can't return history
  control <- pmcmc_control(10, rerun_every = 2, save_trajectories = FALSE,
                           save_state = FALSE)
  ans <- pmcmc(dat$pars, dat$filter, control = control)

  mockery::expect_called(dat$filter$run, 16)
})


test_that("can partially run the pmcmc", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  pars <- pmcmc_parameters$new(
    list(pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                         prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal = proposal_kernel)

  control <- pmcmc_control(30, save_state = TRUE, save_trajectories = TRUE)

  dat <- example_sir()
  n_particles <- 42
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)

  set.seed(1)
  results1 <- pmcmc(pars, p1, control = control)

  control$n_steps_each <- 10
  set.seed(1)
  initial <- pmcmc_check_initial(NULL, pars, 1)[, 1]
  obj <- pmcmc_state$new(pars, initial, p2, control)
  expect_equal(obj$run(), list(step = 10, finished = FALSE))
  tmp <- r6_private(obj)$history_pars$get()
  expect_equal(lengths(tmp), rep(c(2, 0), c(10, 20)))
  expect_equal(obj$run(), list(step = 20, finished = FALSE))
  expect_equal(obj$run(), list(step = 30, finished = TRUE))
  expect_equal(obj$run(), list(step = 30, finished = TRUE))
  results2 <- obj$finish()

  expect_equal(results2, results1)
})


test_that("can change the number of threads mid-run", {
  skip_on_cran()
  dat <- example_sir()
  n_particles <- 20
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)


  control <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE)

  set.seed(1)
  results1 <- pmcmc(dat$pars, p1, control = control)

  control$n_steps_each <- 10
  set.seed(1)
  initial <- pmcmc_check_initial(NULL, dat$pars, 1)[, 1]
  obj <- pmcmc_state$new(dat$pars, initial, p2, control)
  expect_equal(obj$run(), list(step = 10, finished = FALSE))
  expect_equal(obj$set_n_threads(2), 1)
  expect_equal(obj$run(), list(step = 20, finished = FALSE))
  expect_equal(r6_private(r6_private(obj)$filter)$n_threads, 2)
  expect_equal(obj$set_n_threads(1), 2)
  expect_equal(obj$run(), list(step = 30, finished = TRUE))
  expect_equal(r6_private(r6_private(obj)$filter)$n_threads, 1)
  expect_equal(obj$run(), list(step = 30, finished = TRUE))
  results2 <- obj$finish()

  expect_equal(results2, results1)
})


test_that("Can override thread count via control", {
  dat <- example_sir()
  n_particles <- 30
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L, n_threads = 4)
  res <- pmcmc(dat$pars, p,
               control = pmcmc_control(5, n_threads_total = 2))
  expect_equal(res$predict$filter$n_threads, 2)

  p$set_n_threads(4)
  res <- pmcmc(dat$pars, p,
               control = pmcmc_control(5, n_threads_total = NULL))
  expect_equal(res$predict$filter$n_threads, 4)
})


test_that("Can save intermediate state to restart", {
  dat <- example_sir()
  n_particles <- 42
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p3 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  control1 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE)
  control2 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE,
                            save_restart = 20)
  control3 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE,
                            save_restart = c(20, 30))
  set.seed(1)
  res1 <- pmcmc(dat$pars, p1, control = control1)
  set.seed(1)
  res2 <- pmcmc(dat$pars, p2, control = control2)
  set.seed(1)
  res3 <- pmcmc(dat$pars, p3, control = control3)

  ## Same actual run
  expect_identical(res1$trajectories, res2$trajectories)
  expect_identical(res1$trajectories, res3$trajectories)

  expect_null(res1$restart)

  expect_is(res2$restart, "list")
  expect_equal(res2$restart$time, 20)
  expect_equal(dim(res2$restart$state), c(5, 30, 1))

  expect_is(res3$restart, "list")
  expect_equal(res3$restart$time, c(20, 30))
  expect_equal(dim(res3$restart$state), c(5, 30, 2))

  expect_equal(res3$restart$state[, , 1], res2$restart$state[, , 1])
})


test_that("can restart the mcmc using saved state", {
  ## This is mostly a test of "can we" start the pmcmc again, rather
  ## than "should we", and "does it mean anything statistically" as
  ## these are harder questions.
  dat <- example_sir()
  n_particles <- 100

  set.seed(1)
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  control1 <- pmcmc_control(50, save_trajectories = TRUE, save_state = TRUE,
                            progress = FALSE, save_restart = 40)
  res1 <- pmcmc(dat$pars, p1, control = control1)

  ## Our new restart state, which includes a range of possible S
  ## values
  expect_equal(dim(res1$restart$state), c(5, 50, 1))
  s <- res1$restart$state[, , 1]
  d2 <- dat$data[dat$data$day_start >= 40, ]

  initial2 <- particle_filter_initial(s)
  p2 <- particle_filter$new(d2, dat$model, n_particles, dat$compare,
                            index = dat$index, initial = initial2)
  control2 <- pmcmc_control(50, save_trajectories = TRUE, save_state = TRUE,
                            progress = FALSE)
  res2 <- pmcmc(dat$pars, p2, control = control2)

  expect_equal(res2$trajectories$step, (40:100) * 4)
  expect_equal(dim(res2$trajectories$state), c(3, 50, 61))
})


test_that("Fix parameters in sir model", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  pars <- pmcmc_parameters$new(
    list(pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                         prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal = proposal_kernel)

  pars2 <- pars$fix(c(gamma = 0.1))

  dat <- example_sir()
  n_particles <- 40
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  control <- pmcmc_control(10, save_trajectories = TRUE, save_state = TRUE)

  results <- pmcmc(pars2, p, control = control)
  expect_equal(dim(results$pars), c(10, 1))
  expect_equal(results$predict$transform(pi), list(beta = pi, gamma = 0.1))
})


test_that("Split chain manually", {
  dat <- example_sir()
  n_particles <- 30
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  control <- pmcmc_control(10, n_chains = 4, use_parallel_seed = TRUE)

  res1 <- pmcmc(dat$pars, p1, control = control)

  inputs <- pmcmc_chains_prepare(dat$pars, p2, NULL, control)
  samples <- lapply(seq_len(control$n_chains), pmcmc_chains_run, inputs)
  res2 <- pmcmc_combine(samples = samples)

  expect_equal(res1, res2)
})


test_that("split chain running requires one worker", {
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


test_that("Can run pmcmc with early exit enabled", {
  dat <- example_sir()
  ## Really low particle number to inflate variance and cause more
  ## serious drop-outs
  n_particles <- 5
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  control1 <- pmcmc_control(30)
  control2 <- pmcmc_control(30, filter_early_exit = TRUE)
  set.seed(1)
  res1 <- pmcmc(dat$pars, p1, control = control1)
  set.seed(1)
  res2 <- pmcmc(dat$pars, p2, control = control2)

  ## Really hard to test this, the final parameters turn out to be
  ## identical to 300 steps at least, which is surprising. However,
  ## you can see that the RNG streams differ:
  expect_false(identical(p2$inputs()$seed, p1$inputs()$seed))
})


test_that("Fix impossible control parameters", {
  ctrl <- pmcmc_control(10, n_workers = 1)
  ctrl$n_steps_each <- 2
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, 20, dat$compare,
                           index = dat$index, seed = 1L)
  ## Previously this errored, here we're just looking for completion
  results <- pmcmc(dat$pars, p, control = ctrl)
  expect_equal(dim(results$pars), c(10, 2))
  expect_false(any(is.na(results$pars)))
})


test_that("Can save intermediate state to restart", {
  dat <- example_sir()
  p1 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                   index = dat$index)
  p2 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                   index = dat$index)
  p3 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                   index = dat$index)
  control1 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE)
  control2 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE,
                            save_restart = 20)
  control3 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE,
                            save_restart = c(20, 30))

  set.seed(1)
  res1 <- pmcmc(dat$pars, p1, control = control1)
  set.seed(1)
  res2 <- pmcmc(dat$pars, p2, control = control2)
  set.seed(1)
  res3 <- pmcmc(dat$pars, p3, control = control3)

  ## Same actual run
  expect_identical(res1$trajectories, res2$trajectories)
  expect_identical(res1$trajectories, res3$trajectories)

  expect_null(res1$restart)

  expect_is(res2$restart, "list")
  expect_equal(res2$restart$time, 20)
  expect_equal(dim(res2$restart$state), c(5, 30, 1))

  expect_is(res3$restart, "list")
  expect_equal(res3$restart$time, c(20, 30))
  expect_equal(dim(res3$restart$state), c(5, 30, 2))

  expect_equal(res3$restart$state[, , 1], res2$restart$state[, , 1])
})


test_that("Can create restart initial function", {
  m <- matrix(1:70 - 1, 10)
  f <- particle_filter_initial(m)
  expect_identical(parent.env(environment(f)), baseenv())
  res <- f(NULL, 3, NULL)
  expect_equal(dim(res), c(10, 3))
  expect_equal(res %% 10, matrix(0:9, 10, 3))
  expect_true(all(diff(res %/% 10) == 0))
})


test_that("Can filter pmcmc on creation", {
  dat <- example_sir()
  n_particles <- 10
  control1 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE)
  control2 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE,
                            n_burnin = 5, n_steps_retain = 7)

  set.seed(1)
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  results1 <- pmcmc(dat$pars, p, control = control1)
  set.seed(1)
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  results2 <- pmcmc(dat$pars, p, control = control2)

  expect_equal(dim(results1$pars), c(30, 2))
  expect_equal(dim(results2$pars), c(7, 2))

  expect_identical(results1$pars_full, results1$pars)
  expect_identical(results1$probabilities_full, results1$probabilities)
  expect_equal(results2$pars_full, results1$pars)
  expect_equal(results2$probabilities_full, results1$probabilities)

  expect_equal(results1$iteration, 1:30)
  i <- seq(control2$n_burnin + 1,
           by = control2$n_steps_every,
           length.out = control2$n_steps_retain)
  expect_equal(results2$iteration, i)

  cmp <- pmcmc_filter(results1, i)
  v <- setdiff(names(results2), c("pars", "probabilities", "pars_index"))
  expect_equal(results2[v], cmp[v])
  expect_equal(results2$pars, unclass(cmp)$pars)
  expect_equal(results2$probabilities, unclass(cmp)$probabilities)
  expect_null(cmp$pars_index)
})


test_that("Can filter pmcmc on creation, after combining chains", {
  dat <- example_sir()

  n_particles <- 10
  control1 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE,
                            n_chains = 2)
  control2 <- pmcmc_control(30, save_trajectories = TRUE, save_state = TRUE,
                            n_chains = 2, n_burnin = 3, n_steps_retain = 7)

  set.seed(1)
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  results1 <- pmcmc(dat$pars, p, control = control1)
  set.seed(1)
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  results2 <- pmcmc(dat$pars, p, control = control2)

  expect_equal(dim(results2$pars), c(14, 2))

  expect_identical(results1$pars_full, results1$pars)
  expect_identical(results1$probabilities_full, results1$probabilities)

  expect_equal(results2$pars_full, results1$pars)
  expect_equal(results2$probabilities_full, results1$probabilities)

  expect_equal(results1$iteration, rep(1:30, 2))
  i <- seq(control2$n_burnin + 1,
           by = control2$n_steps_every,
           length.out = control2$n_steps_retain)
  expect_equal(results2$iteration, rep(i, 2))

  cmp <- pmcmc_filter(results1, c(i, i + control1$n_steps))
  v <- setdiff(names(results2), c("pars", "probabilities", "pars_index"))
  expect_equal(results2[v], cmp[v])
  expect_equal(results2$pars, unclass(cmp)$pars)
  expect_equal(results2$probabilities, unclass(cmp)$probabilities)
  expect_null(cmp$pars_index)

  results3 <- pmcmc_thin(results2)
  expect_identical(results3$pars, results2$pars)
  expect_identical(results3$pars_full, results3$pars_full)
})
