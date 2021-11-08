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


test_that("rmvnorm_generator skips zero variables", {
  n <- 10
  vcv <- cov(matrix(rnorm(n * n), n, n))
  i <- c(2, 5)
  vcv[i, ] <- 0
  vcv[, i] <- 0

  ## mvtnorm does poorly here:
  ## > cmp <- replicate(100, drop(mvtnorm::rmvnorm(1, sigma = vcv)))[i, ]

  set.seed(1)
  f <- rmvnorm_generator(vcv)
  res <- replicate(100, f(rep(0, n)))
  expect_equal(res[i, ], matrix(0, 2, 100))
})


test_that("rmvnorm_generator can created scaled samples", {
  vcv <- matrix(c(4, 2, 2, 3), ncol = 2)
  theta <- runif(2)
  set.seed(1)
  y1 <- rmvnorm_generator(vcv)(theta)
  set.seed(1)
  y2 <- rmvnorm_generator(vcv)(theta, 2)
  set.seed(1)
  y3 <- rmvnorm_generator(vcv * 2)(theta)

  expect_equal((y1 - theta) * sqrt(2), y2 - theta)
  expect_equal(y2, y3)
})

test_that("NLAYER", {
  expect_equal(NLAYER(array(1, c(2, 2, 2))), 2)
  expect_equal(NLAYER(1), 1)
})

test_that("layernames", {
  expect_equal(layernames(array(1, c(2, 2, 2))), NULL)
  expect_equal(layernames(array(1, c(1, 1, 1), as.list(letters[1:3]))), "c")
})

test_that("layernames<-", {
  x <- matrix(1)
  expect_error({
    layernames(x) <- "a"
  },
  "less than")

  x1 <- x2 <- array(1, c(1, 1, 1))
  expect_error({
    layernames(x2) <- NULL
  }, "cannot be")

  x <- array(1, c(1, 1, 1))
  layernames(x) <- "a"
  expect_equal(x, array(1, c(1, 1, 1), list(NULL, NULL, "a")))

  x <- array(1, c(1, 1, 1), dimnames = list(NULL, NULL, "b"))
  layernames(x) <- "a"
  expect_equal(x, array(1, c(1, 1, 1), list(NULL, NULL, "a")))

  x <- array(1, c(1, 1, 1), dimnames = list(NULL, NULL, "b"))
  layernames(x) <- NULL
  expect_equal(x, array(1, c(1, 1, 1), list(NULL, NULL, NULL)))
})

test_that("set_layernames", {
  x <- (array(1, c(1, 1, 1), as.list(letters[1:3])))
  expect_equal(set_layernames(array(1, c(1, 1, 1), list("a", "b", NULL)), "c"),
               x)
})

test_that("is_3d_array", {
  expect_true(is_3d_array(array(1, c(1, 1, 1))))
  expect_false(is_3d_array(matrix(1, 2, 2)))
})
