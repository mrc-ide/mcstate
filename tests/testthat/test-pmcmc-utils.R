context("pmcmc (utils)")

test_that("format and print the simplest object", {
  pars <- matrix(NA_real_, 10, 4,
                 dimnames = list(NULL, c("a", "b", "c", "d")))
  probs <- matrix(NA_real_, 10, 3, dimnames = list(NULL, c("x", "y", "z")))
  x <- mcstate_pmcmc(pars, probs, NULL, NULL, NULL, NULL)

  expected <- c(
    "<mcstate_pmcmc> (10 samples)",
    "  pars: 10 x 4 matrix of parameters",
    "    a, b, c, d",
    "  probabilities: 10 x 3 matrix of log-probabilities",
    "    x, y, z",
    "  state: (not included)",
    "  trajectories: (not included)",
    "  restart: (not included)")

  expect_equal(format(x), expected)
  expect_output(print(x), paste(expected, collapse = "\n"), fixed = TRUE)
})


test_that("format and print with state", {
  pars <- matrix(NA_real_, 10, 4,
                 dimnames = list(NULL, c("a", "b", "c", "d")))
  probs <- matrix(NA_real_, 10, 3, dimnames = list(NULL, c("x", "y", "z")))
  state <- matrix(NA_real_, 4, 10)
  trajectories <- list(state = array(NA_real_, c(4, 10, 20)))
  predict <- NULL
  restart <- list(date = 1, state = array(NA_real_, c(4, 10, 1)))

  x <- mcstate_pmcmc(pars, probs, state, trajectories, restart, predict)

  expected <- c(
    "<mcstate_pmcmc> (10 samples)",
    "  pars: 10 x 4 matrix of parameters",
    "    a, b, c, d",
    "  probabilities: 10 x 3 matrix of log-probabilities",
    "    x, y, z",
    "  state: 4 x 10 matrix of final states",
    "  trajectories: 4 x 10 x 20 array of particle trajectories",
    "  restart: 4 x 10 x 1 array of particle restart state")

  expect_equal(format(x), expected)
  expect_output(print(x), paste(expected, collapse = "\n"), fixed = TRUE)
})


test_that("print multichain object", {
  x <- pmcmc_combine(samples = example_sir_pmcmc2()$results)

  expected <- c(
    "<mcstate_pmcmc> (93 samples across 3 chains)",
    "  pars: 93 x 2 matrix of parameters",
    "    beta, gamma",
    "  probabilities: 93 x 3 matrix of log-probabilities",
    "    log_prior, log_likelihood, log_posterior",
    "  state: 5 x 93 matrix of final states",
    "  trajectories: 3 x 93 x 101 array of particle trajectories",
    "  restart: 5 x 93 x 1 array of particle restart state")

  expect_equal(format(x), expected)
  expect_output(print(x), paste(expected, collapse = "\n"), fixed = TRUE)
})


test_that("wrap long variable names nicely", {
  nms <- strrep(c("a", "b", "c", "d"), c(10, 20, 30, 40))
  pars <- matrix(NA_real_, 10, 4, dimnames = list(NULL, nms))
  probs <- matrix(NA_real_, 10, 3, dimnames = list(NULL, c("x", "y", "z")))
  x <- withr::with_options(
    list(width = 80),
    format(mcstate_pmcmc(pars, probs, NULL, NULL, NULL, NULL)))
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


test_that("progress bar is a noop when progress = FALSE", {
  p <- pmcmc_progress(3, FALSE, force = TRUE)
  expect_silent(p())
  expect_null(p())
  expect_silent(p())
})


test_that("progress bar creates progress_bar when progress = TRUE", {
  p <- pmcmc_progress(3, TRUE, force = TRUE)
  Sys.sleep(0.2)
  expect_message(
    p(),
    "Step 1 / 3 \\[=*>-+\\] ETA .* \\| 00:00:[0-9]{2} so far")
  expect_s3_class(
    suppressMessages(p()),
    "progress_bar")
  expect_message(
    p(),
    "Finished 3 steps in [0-9]+ secs")
})
