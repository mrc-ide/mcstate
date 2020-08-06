context("pmcmc_parameters")

test_that("Can construct a parameter", {
  p <- pmcmc_parameter(1, 0, 10)
  expect_s3_class(p, "pmcmc_parameter")
  expect_equal(p$initial, 1)
  expect_equal(p$min, 0)
  expect_equal(p$max, 10)
  expect_false(p$discrete)
  expect_equal(p$prior(1), 0)
})


test_that("Can provide a prior", {
  f <- function(p) log(1/p)
  p <- pmcmc_parameter(1, prior = f)
  expect_identical(p$prior, f)
})


## TODO: Break this up a bit
test_that("parameters", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  p <- pmcmc_parameters$new(
    list(
      beta = pmcmc_parameter(0.2, min = 0, max = 1,
                             prior = function(p) log(1e-10)),
      gamma = pmcmc_parameter(0.1, min = 0, max = 1,
                              prior = function(p) log(1e-10))),
    proposal = proposal_kernel)

  expect_equal(
    p$initial(),
    c(beta = 0.2, gamma = 0.1))

  set.seed(1)
  p0 <- p$initial()
  p1 <- p$propose(p0)
  set.seed(1)
  expect_equal(p1, drop(rmvnorm(1, p0, proposal_kernel)))
})


test_that("transform parameters returns a list by default", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  p <- pmcmc_parameters$new(
    list(
      beta = pmcmc_parameter(0.2, min = 0, max = 1,
                             prior = function(p) log(1e-10)),
      gamma = pmcmc_parameter(0.1, min = 0, max = 1,
                              prior = function(p) log(1e-10))),
    proposal = proposal_kernel)

  expect_equal(
    p$model(p$initial()),
    list(beta = 0.2, gamma = 0.1))
  expect_equal(
    p$model(p$initial() + 0.5),
    list(beta = 0.7, gamma = 0.6))
})


test_that("provide custom transform and use it", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  transform <- function(p) {
    list(alpha = p[["beta"]] + p[["gamma"]],
         gamma = rep(p[["gamma"]], 4),
         beta = p[["beta"]])
  }

  p <- pmcmc_parameters$new(
    list(beta = pmcmc_parameter(0.2), gamma = pmcmc_parameter(0.1)),
    proposal = proposal_kernel,
    transform = transform)

  expect_equal(
    p$model(p$initial()),
    list(alpha = 0.3, gamma = rep(0.1, 4), beta = 0.2))
  expect_equal(
    p$model(p$initial() + 0.5),
    list(alpha = 1.3, gamma = rep(0.6, 4), beta = 0.7))
})


test_that("can compute prior", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  transform <- function(p) {
    list(alpha = p[["beta"]] + p[["gamma"]],
         gamma = rep(p[["gamma"]], 4),
         beta = p[["beta"]])
  }

  p <- pmcmc_parameters$new(
    list(
      beta = pmcmc_parameter(0.2, prior = function(p) dexp(p, log = TRUE)),
      gamma = pmcmc_parameter(0.1, prior = function(p) dnorm(p, log = TRUE))),
    proposal = diag(2),
    transform = transform)

  ## Compute prior
  expect_equal(
    p$prior(p$initial()),
    dexp(0.2, log = TRUE) + dnorm(0.1, log = TRUE))
  ## Done by position, not name:
  expect_equal(
    p$prior(list(0.2, 0.1)),
    dexp(0.2, log = TRUE) + dnorm(0.1, log = TRUE))
  ## Non-default values
  expect_equal(
    p$prior(list(beta = 1, gamma = 0)),
    dexp(1, log = TRUE) + dnorm(0, log = TRUE))
})
