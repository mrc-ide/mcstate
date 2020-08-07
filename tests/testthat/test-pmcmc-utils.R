context("pmcmc (utils)")

test_that("format and print the simplest object", {
  pars <- matrix(NA_real_, 10, 4,
                 dimnames = list(NULL, c("a", "b", "c", "d")))
  probs <- matrix(NA_real_, 10, 3, dimnames = list(NULL, c("x", "y", "z")))
  x <- mcstate_pmcmc(pars, probs, NULL, NULL)

  expected <- c(
    "<mcstate_pmcmc> (9 samples)",
    "  pars: 10 x 4 matrix of parameters",
    "    a, b, c, d",
    "  probabilities: 10 x 3 matrix of log-probabilities",
    "    x, y, z",
    "  state: (not included)",
    "  trajectories: (not included)")

  expect_equal(format(x), expected)
  expect_output(print(x), paste(expected, collapse = "\n"), fixed = TRUE)
})


test_that("format and print with state", {
  pars <- matrix(NA_real_, 10, 4,
                 dimnames = list(NULL, c("a", "b", "c", "d")))
  probs <- matrix(NA_real_, 10, 3, dimnames = list(NULL, c("x", "y", "z")))
  state <- matrix(NA_real_, 4, 10)
  trajectories <- array(NA_real_, c(4, 20, 10))

  x <- mcstate_pmcmc(pars, probs, state, trajectories)

  expected <- c(
    "<mcstate_pmcmc> (9 samples)",
    "  pars: 10 x 4 matrix of parameters",
    "    a, b, c, d",
    "  probabilities: 10 x 3 matrix of log-probabilities",
    "    x, y, z",
    "  state: 4 x 10 matrix of final states",
    "  trajectories: 4 x 20 x 10 array of particle trajectories")

  expect_equal(format(x), expected)
  expect_output(print(x), paste(expected, collapse = "\n"), fixed = TRUE)
})


test_that("wrap long variable names nicely", {
  nms <- strrep(c("a", "b", "c", "d"), c(10, 20, 30, 40))
  pars <- matrix(NA_real_, 10, 4, dimnames = list(NULL, nms))
  probs <- matrix(NA_real_, 10, 3, dimnames = list(NULL, c("x", "y", "z")))
  x <- withr::with_options(
    list(width = 80),
    format(mcstate_pmcmc(pars, probs, NULL, NULL)))
  expect_equal(
    x[[3]],
    "    aaaaaaaaaa, bbbbbbbbbbbbbbbbbbbb, cccccccccccccccccccccccccccccc,")
  expect_equal(
    x[[4]],
    "    dddddddddddddddddddddddddddddddddddddddd")
})


test_that("can filter trajectories with dropping dimensions", {
  m <- array(1:20, c(1, 4, 5))
  expect_equal(
    sample_trajectory(m, 2),
    matrix(m[, 2, ], 1, 5))
  expect_equal(
    sample_trajectory(m, 2:3),
    matrix(m[, 2:3, ], 2, 5))
})
