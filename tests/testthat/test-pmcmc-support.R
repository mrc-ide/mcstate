context("pmcmc (support)")

test_that("rerun determinstically", {
  i <- 0:10
  every <- 3
  set.seed(1)
  x <- replicate(1000, vlapply(i, rerun, 3, FALSE))
  y <- replicate(1000, vlapply(i, rerun, 3, TRUE))


  expect_equal(x, matrix(rep(c(TRUE, FALSE, FALSE), length.out = 11),
                         11, 1000))
  expect_false(all(x == y))
  expect_gt(mean(y), 0.32)
  expect_lt(mean(y), 0.34)

  expect_equal(replicate(10, vlapply(i, rerun, Inf, FALSE)),
               matrix(FALSE, 11, 10))
  expect_equal(replicate(10, vlapply(i, rerun, Inf, TRUE)),
               matrix(FALSE, 11, 10))
})
