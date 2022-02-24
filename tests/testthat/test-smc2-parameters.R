context("smc2_parameters")

## This is where the issue is at the moment. What we really want to do
## is accept a d and an r function here, along with parameters for
## them. For the sampling we need to accept a number of samples
test_that("can construct a parameter", {
  p <- smc2_parameter("a",
                     function(n) runif(n, 0, 10),
                     function(x) dunif(x, 0, 10, log = TRUE),
                     min = 0, max = 10, integer = TRUE)

  expect_s3_class(p, "smc2_parameter")
  expect_equal(p$name, "a")
  expect_equal(p$min, 0)
  expect_equal(p$max, 10)
  expect_true(p$integer)
  expect_equal(p$prior(1), log(0.1))
  set.seed(1)
  res <- p$sample(20)
  set.seed(1)
  expect_equal(res, runif(20, 0, 10))
})


test_that("Can use 'discrete' argument but deprecation warning is shown", {
  expect_warning(p <- smc2_parameter("a",
                                     function(n) runif(n, 0, 10),
                                     function(x) dunif(x, 0, 10, log = TRUE),
                                     min = 0, max = 10, discrete = TRUE),
                 "'discrete' is deprecated.\nUse 'integer' instead.")
  expect_s3_class(p, "smc2_parameter")
  expect_equal(p$integer, TRUE)
})


test_that("smc2_parmeter must satisfy min/max constraints", {
  expect_error(
    smc2_parameter("a",
                   function(n) runif(n, 0, 10),
                   function(x) dunif(x, 0, 10, log = TRUE),
                   min = 20, max = 10),
    "'max' must be > 'min' (20)",
    fixed = TRUE)
})


test_that("can construct a set of parameters", {
  p <- smc2_parameters$new(
    list(
      smc2_parameter("beta",
                     function(n) runif(n, 0, 10),
                     function(x) dunif(x, 0, 10, log = TRUE),
                     min = 0, max = 10),
      smc2_parameter("gamma",
                     function(n) runif(n, 1, 2),
                     function(x) dunif(x, 1, 2, log = TRUE),
                     min = 1, max = 2)))
  expect_s3_class(p, "smc2_parameters")
  expect_equal(p$names(), c("beta", "gamma"))
  expect_equal(p$summary(),
               data_frame(name = c("beta", "gamma"),
                          min = c(0, 1),
                          max = c(10, 2),
                          discrete = FALSE,
                          integer = FALSE))
})


test_that("can sample initial conditions from the prior", {
  p <- smc2_parameters$new(
    list(
      smc2_parameter("beta",
                     function(n) runif(n, 0, 10),
                     function(x) dunif(x, 0, 10, log = TRUE),
                     min = 0, max = 10),
      smc2_parameter("gamma",
                     function(n) runif(n, 1, 2),
                     function(x) dunif(x, 1, 2, log = TRUE),
                     min = 1, max = 2)))

  theta <- p$sample(10)
  ## expect_is deprecated, but no other way of doing this directly
  ## with testthat now: https://github.com/r-lib/testthat/issues/865
  expect_true(is.matrix(theta))
  expect_equal(colnames(theta), c("beta", "gamma"))
  expect_equal(dim(theta), c(10, 2))
})


test_that("can propose new parameters", {
  p <- smc2_parameters$new(
    list(
      smc2_parameter("beta",
                     function(n) runif(n, 0, 10),
                     function(x) dunif(x, 0, 10, log = TRUE),
                     min = 0, max = 10),
      smc2_parameter("gamma",
                     function(n) rexp(n),
                     function(x) dexp(x, log = TRUE))))

  theta <- p$sample(100)
  vcv <- diag(2) * 100

  new <- p$propose(theta, vcv)

  expect_true(is.matrix(new))
  expect_equal(colnames(new), c("beta", "gamma"))
  expect_equal(dim(new), c(100, 2))
  expect_true(all(new[, 1] > 0 & new[, 1] < 10))
})
