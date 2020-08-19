context("pmcmc_parameters")

test_that("Can construct a parameter", {
  p <- pmcmc_parameter("a", 1, 0, 10)
  expect_s3_class(p, "pmcmc_parameter")
  expect_equal(p$name, "a")
  expect_equal(p$initial, 1)
  expect_equal(p$min, 0)
  expect_equal(p$max, 10)
  expect_false(p$discrete)
  expect_equal(p$prior(1), 0)
})


test_that("Can provide a prior", {
  f <- function(p) log(1 / p)
  p <- pmcmc_parameter("a", 1, prior = f)
  expect_identical(p$prior, f)
})


test_that("parameter initial point must lie in range", {
  expect_silent(pmcmc_parameter("a", 0))
  expect_silent(pmcmc_parameter("a", 0, min = 0))
  expect_silent(pmcmc_parameter("a", 0, max = 0))
  expect_error(pmcmc_parameter("a", 0, min = 1),
               "'initial' must be >= 'min' (1)", fixed = TRUE)
  expect_error(pmcmc_parameter("a", 0, max = -1),
               "'initial' must be <= 'max' (-1)", fixed = TRUE)
})


test_that("initial value must satify prior and not fail", {
  expect_silent(
    pmcmc_parameter("a", 10, prior = log))
  expect_error(
    pmcmc_parameter("a", 0, prior = log),
    "Prior function for 'a' returned a non-finite value on initial value")
  expect_error(
    pmcmc_parameter("a", 0, prior = function(p) stop("an error")),
    "Prior function for 'a' failed to evaluate initial value: an error")
})


test_that("parameters", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  p <- pmcmc_parameters$new(
    list(
      pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                      prior = function(p) log(1e-10)),
      pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                      prior = function(p) log(1e-10))),
    proposal = proposal_kernel)

  expect_equal(
    p$initial(),
    c(beta = 0.2, gamma = 0.1))

  set.seed(1)
  p0 <- p$initial()
  p1 <- p$propose(p0)
  set.seed(1)
  expect_equal(p1, rmvnorm_generator(proposal_kernel)(p0))
})


test_that("transform parameters returns a list by default", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  p <- pmcmc_parameters$new(
    list(
      pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                      prior = function(p) log(1e-10)),
      pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
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
    list(pmcmc_parameter("beta", 0.2), pmcmc_parameter("gamma", 0.1)),
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
      pmcmc_parameter("beta", 0.2, prior = function(p) dexp(p, log = TRUE)),
      pmcmc_parameter("gamma", 0.1, prior = function(p) dnorm(p, log = TRUE))),
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


test_that("all inputs to pmcmc_parameters must be pmcmc_parameter objects", {
  expect_error(
    pmcmc_parameters$new(NULL, proposal = diag(2)),
    "'parameters' must be a list")
  expect_error(
    pmcmc_parameters$new(list(), proposal = diag(2)),
    "At least one parameter is required")

  expect_error(
    pmcmc_parameters$new(list(a = 1, b = 2), proposal = diag(2)),
    "Expected all elements of '...' to be 'pmcmc_parameter' objects",
    fixed = TRUE)
})


test_that("proposal matrix and parameters must conform", {
  pars <- list(pmcmc_parameter("beta", 0.2), pmcmc_parameter("gamma", 0.1))
  expect_error(
    pmcmc_parameters$new(pars, proposal = diag(1)),
    "Expected a square proposal matrix with 2 rows and columns")
  expect_error(
    pmcmc_parameters$new(pars, proposal = diag(3)),
    "Expected a square proposal matrix with 2 rows and columns")
  expect_error(
    pmcmc_parameters$new(pars, proposal = matrix(0, 3, 2)),
    "Expected a square proposal matrix with 2 rows and columns")
})


test_that("if proposal has names, they must match parameters", {
  pars <- list(pmcmc_parameter("beta", 0.2), pmcmc_parameter("gamma", 0.1))
  nms <- c("beta", "gamma")
  m <- function(a, b = a) {
    matrix(1, length(a), length(b), dimnames = list(a, b))
  }

  expect_error(
    pmcmc_parameters$new(pars, proposal = m(c("a", "b"))),
    "Expected dimension names of 'proposal' to match parmeters")
  expect_error(
    pmcmc_parameters$new(pars, proposal = m(c("a", "b"), nms)),
    "Expected dimension names of 'proposal' to match parmeters")
  expect_error(
    pmcmc_parameters$new(pars, proposal = m(nms, c("a", "b"))),
    "Expected dimension names of 'proposal' to match parmeters")
  expect_silent(pmcmc_parameters$new(pars, proposal = m(nms)))
})


test_that("can summarise parameters", {
  p <- pmcmc_parameters$new(
    list(
      pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                      prior = function(p) log(1e-10)),
      pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                      prior = function(p) log(1e-10))),
    proposal = diag(2))

  expect_equal(p$names(), c("beta", "gamma"))
  expect_equal(
    p$summary(),
    data_frame(name = c("beta", "gamma"), min = 0, max = 1, discrete = FALSE))
})
