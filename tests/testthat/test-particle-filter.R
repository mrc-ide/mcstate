context("particle_filter")

test_that("run particle filter on sir model", {
  dat <- example_sir()

  p <- particle_filter$new(dat$data, dat$compare, FALSE)
  res <- p$run(dat$y0, dat$model, 42)
  expect_is(res, "numeric")

  expect_is(p$state, "matrix")
  expect_equal(dim(p$state), c(3, 42))
  expect_null(p$history)
})


test_that("particle filter likelihood is worse with worse parameters", {
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$compare, FALSE)
  ll1 <- p$run(dat$y0, dat$model(), 100)
  ll2 <- p$run(dat$y0, dat$model(gamma = 1, beta = 1), 100)
  expect_true(ll1 > ll2)
})
