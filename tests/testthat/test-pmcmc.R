context("pmcmc")


## TODO: Ed and Lilith we could use better tests throughout here. The
## sampler should be ok to run for ~10k iterations without taking too
## long to be annoying in tests.
test_that("mcmc works for uniform distribution on unit square", {
  dat <- example_uniform()

  set.seed(1)
  testthat::try_again(5, {
    res <- pmcmc(dat$pars, dat$filter, 1000, FALSE, FALSE)
    expect_s3_class(res, "mcstate_pmcmc")
    expect_true(all(acceptance_rate(res$pars) == 1))
    expect_true(abs(mean(res$pars[, "a"]) - 0.5) < 0.05)
    expect_true(abs(mean(res$pars[, "b"]) - 0.5) < 0.05)
  })
})


test_that("mcmc works for multivariate gaussian", {
  dat <- example_mvnorm()

  set.seed(1)
  testthat::try_again(5, {
    res <- pmcmc(dat$pars, dat$filter, 1000, FALSE, FALSE)
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
  res <- pmcmc(dat$pars, dat$filter, 100, FALSE, FALSE)

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

  set.seed(1)
  results1 <- pmcmc(pars, p1, 30, TRUE, TRUE)
  set.seed(1)
  results2 <- pmcmc(pars, p2, 30, FALSE, FALSE)

  expect_setequal(
    names(results1),
    c("chain", "iteration",
      "pars", "probabilities", "state", "trajectories", "predict"))

  expect_null(results1$chain)
  expect_equal(results1$iteration, 0:30)

  ## Including or not the history does not change the mcmc trajectory:
  expect_identical(names(results1), names(results2))
  expect_equal(results1$pars, results2$pars)
  expect_equal(results1$probabilities, results2$probabilities)

  ## We did mix
  expect_true(all(acceptance_rate(results1$pars) > 0))

  ## Parameters and probabilities have the expected shape
  expect_equal(dim(results1$pars), c(31, 2))
  expect_equal(colnames(results1$pars), c("beta", "gamma"))

  expect_equal(dim(results1$probabilities), c(31, 3))
  expect_equal(colnames(results1$probabilities),
               c("log_prior", "log_likelihood", "log_posterior"))

  ## History, if returned, has the correct shape
  expect_equal(dim(results1$state), c(4, 31)) # state, mcmc

  ## Trajectories, if returned, have the same shape
  expect_s3_class(results1$trajectories, "mcstate_trajectories")
  expect_equal(dim(results1$trajectories$state), c(3, 31, 101))
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
  results1 <- pmcmc(dat$pars, p1, 5, FALSE, FALSE)
  set.seed(1)
  results2 <- pmcmc(dat$pars, p2, 5, TRUE, FALSE)
  set.seed(1)
  results3 <- pmcmc(dat$pars, p3, 5, TRUE, FALSE)

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
  expect_message(
    pmcmc(dat$pars, dat$filter, 1000, FALSE, FALSE, progress = TRUE),
    "Finished 1000 steps in ")
})


test_that("run multiple chains", {
  dat <- example_uniform()

  set.seed(1)
  res1 <- pmcmc(dat$pars, dat$filter, 100, FALSE, FALSE, n_chains = 1)
  expect_s3_class(res1, "mcstate_pmcmc")
  expect_null(res1$chain)

  set.seed(1)
  res3 <- pmcmc(dat$pars, dat$filter, 100, FALSE, FALSE, n_chains = 3)
  expect_s3_class(res3, "mcstate_pmcmc")
  expect_equal(res3$chain, rep(1:3, each = 101))

  expect_equal(res1$pars, res3$pars[1:101, ])
})


test_that("progress in multiple chains", {
  dat <- example_uniform()
  expect_message(
    pmcmc(dat$pars, dat$filter, 100, FALSE, FALSE, progress = TRUE,
          n_chains = 3),
    "Running chain 2 / 3")
})


test_that("All arguments forwarded to multiple chains", {
  skip_if_not_installed("mockery")
  dat <- example_uniform()
  res <- list(
    pmcmc(dat$pars, dat$filter, 1000, FALSE, FALSE),
    pmcmc(dat$pars, dat$filter, 1000, FALSE, FALSE),
    pmcmc(dat$pars, dat$filter, 1000, FALSE, FALSE))

  mock_pmcmc_single_chain <- mockery::mock(
    res[[1]], res[[1]], res[[2]], res[[3]])
  with_mock("mcstate::pmcmc_single_chain" = mock_pmcmc_single_chain, {
    ans1 <- pmcmc(dat$pars, dat$filter, 1000, FALSE, FALSE, FALSE)
    ans2 <- pmcmc(dat$pars, dat$filter, 1000, FALSE, FALSE, FALSE, n_chains = 3)
  })

  mockery::expect_called(mock_pmcmc_single_chain, 4)
  args <- mockery::mock_args(mock_pmcmc_single_chain)
  ic <- c(a = 0.5, b = 0.5)
  expect_equal(args[[1]],
               list(dat$pars, ic, dat$filter, 1000, Inf, FALSE, FALSE, FALSE))
  expect_equal(args[[2]],
               list(dat$pars, ic, dat$filter, 1000, Inf, FALSE, FALSE, FALSE))
  expect_equal(args[[3]],
               list(dat$pars, ic, dat$filter, 1000, Inf, FALSE, FALSE, FALSE))
  expect_equal(args[[4]],
               list(dat$pars, ic, dat$filter, 1000, Inf, FALSE, FALSE, FALSE))

  expect_equal(ans1, res[[1]])
  expect_equal(ans2, pmcmc_combine(samples = res))
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

  set.seed(1)
  results1 <- pmcmc(pars, p1, 5, TRUE, TRUE)
  set.seed(1)
  results2 <- pmcmc(pars, p2, 5, TRUE, TRUE)

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
               "Expected a vector of length 2 for 'initial'")
  expect_error(pmcmc_check_initial(c(x = 0.1, y = 0.2), dat$pars, 1),
               "If 'initial' has names, they must match pars$names()",
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
    "Expected a matrix with 2 rows for 'initial'")
  expect_error(
    pmcmc_check_initial(matrix(0.5, 2, 6), dat$pars, 5),
    "Expected a matrix with 5 columns for 'initial'")

  expect_error(
    pmcmc_check_initial(matrix(0.5, 2, 5, dimnames = list(c("x", "y"), NULL)),
                        dat$pars, 5),
    "If 'initial' has rownames, they must match pars$names()",
    fixed = TRUE)

  m <- matrix(runif(10), 2, 5)
  i <- cbind(c(2, 1, 2), c(2, 4, 5))
  m[i] <- -m[i]
  expect_error(
    pmcmc_check_initial(m, dat$pars, 5),
    "Starting point does not have finite prior probability (2, 4, 5)",
    fixed = TRUE)
})


test_that("can start a pmcmc from a matrix of starting points", {
  dat <- example_uniform()
  initial <- matrix(runif(6), 2, 3, dimnames = list(c("a", "b"), NULL))
  res <- pmcmc(dat$pars, dat$filter, 1000, FALSE, FALSE,
               n_chains = 3, initial = initial)
  expect_equal(res$pars[res$iteration == 0, ], t(initial))
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

  set.seed(1)
  res1 <- pmcmc(pars, p1, 30, TRUE, TRUE, rerun_every = 2)
  expect_lte(max(rle(res1$probabilities[, "log_likelihood"])$lengths), 2)

  set.seed(1)
  res2 <- pmcmc(pars, p1, 30, TRUE, TRUE, rerun_every = Inf)
  expect_gt(max(rle(res2$probabilities[, "log_likelihood"])$lengths), 5)
})


test_that("rerunning the particle filter triggers the filter run method", {
  skip_if_not_installed("mockery")
  dat <- example_uniform()
  dat$filter$run <- mockery::mock(1, cycle = TRUE)
  dat$inputs <- function() NULL

  ## with the dummy version we can't return history
  ans <- pmcmc(dat$pars, dat$filter, 10, rerun_every = 2,
               save_trajectories = FALSE, save_state = FALSE)

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

  dat <- example_sir()
  n_particles <- 100
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)

  set.seed(1)
  results1 <- pmcmc(pars, p1, 30, TRUE, TRUE)

  ## This is quite manual here:
  set.seed(1)
  initial <- pmcmc_check_initial(NULL, pars, 1)[, 1]
  obj <- pmcmc_state$new(pars, initial, p2, 30, 10, Inf, TRUE, TRUE, FALSE)
  expect_equal(obj$run(), list(step = 10, finished = FALSE))
  tmp <- r6_private(obj)$history_pars$get()
  expect_equal(lengths(tmp), rep(c(2, 0), c(11, 20)))
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

  set.seed(1)
  results1 <- pmcmc(dat$pars, p1, 30, TRUE, TRUE)

  ## This is quite manual here:
  set.seed(1)
  initial <- pmcmc_check_initial(NULL, dat$pars, 1)[, 1]
  obj <- pmcmc_state$new(dat$pars, initial, p2, 30, 10, Inf, TRUE, TRUE, FALSE)
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
