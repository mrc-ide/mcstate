context("pmcmc (support)")

test_that("rerun determinstically", {
  i <- 0:10
  every <- 3
  set.seed(1)
  rerun_x <- make_rerun(3, FALSE)
  rerun_y <- make_rerun(3, TRUE)

  x <- replicate(1000, vlapply(i, rerun_x))
  y <- replicate(1000, vlapply(i, rerun_y))

  expect_equal(x, matrix(rep(c(TRUE, FALSE, FALSE), length.out = 11),
                         11, 1000))
  expect_false(all(x == y))
  expect_gt(mean(y), 0.32)
  expect_lt(mean(y), 0.34)

  rerun <- make_rerun(Inf, FALSE)
  expect_equal(replicate(10, vlapply(i, rerun)),
                         matrix(FALSE, 11, 10))
  rerun <- make_rerun(Inf, TRUE)
  expect_equal(replicate(10, vlapply(i, rerun)),
               matrix(FALSE, 11, 10))
})
