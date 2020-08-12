context("utils")

test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})


test_that("rmvnorm_generator agrees with rmvnorm", {
  testthat::skip_on_cran() # depends on another package's internals
  vcv <- matrix(c(4, 2, 2, 3), ncol = 2)
  set.seed(1)
  cmp <- replicate(100, drop(mvtnorm::rmvnorm(1, sigma = vcv)))

  set.seed(1)
  f <- rmvnorm_generator(vcv)
  res <- replicate(100, f(c(0, 0)))

  expect_equal(res, cmp)
})


test_that("rmvnorm_generator requires symmetric matrix", {
  expect_error(
    rmvnorm_generator(matrix(1:4, 2, 2)),
    "vcv must be symmetric")
  expect_error(
    rmvnorm_generator(matrix(1, 2, 6)),
    "vcv must be symmetric")
})


test_that("rmvnorm_generator requires positive definite", {
  expect_error(
    rmvnorm_generator(matrix(c(1, 2, 2, 1), 2, 2)),
    "vcv must be positive definite")
})
