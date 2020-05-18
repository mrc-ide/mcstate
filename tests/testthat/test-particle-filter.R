context("particle_filter")

test_that("run particle filter on sir model", {
  dat <- example_sir()

  p <- particle_filter$new(dat$data, dat$compare, FALSE)
  res <- p$run(dat$y0, dat$model(), 42)
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


test_that("can be started with varying initial states", {
  dat <- example_sir()
  n <- sample(10:20, 20, replace = TRUE)
  y0 <- rbind(1010 - n, n, 0, deparse.level = 0)

  ## If we never return a likelihood, the we never shuffle trajectories...
  p <- particle_filter$new(dat$data, function(...) NULL, TRUE)
  ll <- p$run(y0, dat$model(), 20)
  ## ...and therefore we recover our initial states
  expect_equal(p$history[, , 1], y0)
})


test_that("validate particle initial state with vector creates matrix", {
  expect_equal(particle_initial_state(1:5, 10),
               matrix(1:5, 5, 10))
})


test_that("validate particle initial state validates matrix", {
  m <- matrix(runif(12), 3, 4)
  expect_error(particle_initial_state(m, 3),
               "Expected '3' columns for initial state")
  expect_identical(particle_initial_state(m, 4), m)
})


test_that("stop simulation when likelihood is impossible", {
  dat <- example_sir()

  compare <- function(state, output, observed) {
    ret <- dat$compare(state, output, observed)
    if (observed$incid > 15) {
      ret[] <- -Inf
    }
    ret
  }

  p <- particle_filter$new(dat$data, compare, TRUE)
  res <- p$run(dat$y0, dat$model(), 42)
  expect_equal(res, -Inf)

  i <- (which(dat$data$incid > 15)[[1]] + 2):101
  expect_false(any(is.na(p$history[, , !i])))
  expect_true(all(is.na(p$history[, , i])))
})


test_that("Data validation finds correct columns", {
  expect_error(
    particle_filter_validate_data(data.frame(a = 1:10)),
    "Expected columns missing from data: 'step_start', 'step_end'")
  expect_error(
    particle_filter_validate_data(data.frame(step_start = 1:10)),
    "Expected columns missing from data: 'step_end'")
  expect_error(
    particle_filter_validate_data(data.frame(step_end = 1:10)),
    "Expected columns missing from data: 'step_start'")
})


test_that("Data validation requires two times or more", {
  d <- data.frame(step_start = 0, step_end = 10)
  expect_error(
    particle_filter_validate_data(d),
    "Expected at least two time windows")
  expect_error(
    particle_filter_validate_data(d[integer(0), ]),
    "Expected at least two time windows")
})


test_that("Data validation consecutive time windows", {
  d <- data.frame(step_start = c(0, 5, 10),
                  step_end = c(6, 11, 16))
  expect_error(
    particle_filter_validate_data(d),
    "Expected time windows to be adjacent")
})
