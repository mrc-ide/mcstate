context("particle_filter")

test_that("run particle filter on sir model", {
  dat <- example_sir()

  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  res <- p$run(dat$y0, 42, FALSE)
  expect_is(res, "numeric")

  expect_is(p$state, "matrix")
  expect_equal(dim(p$state), c(5, 42))
  expect_null(p$history)
})


test_that("particle filter likelihood is worse with worse parameters", {
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  ll1 <- p$run(dat$y0, 100)
  ll2 <- p$run(dat$y0, 100, user = list(gamma = 1, beta = 1))
  expect_true(ll1 > ll2)
})


test_that("can be started with varying initial states", {
  dat <- example_sir()

  n_particles <- 20
  n <- sample(10:20, n_particles, replace = TRUE)
  y0 <- rbind(1010 - n, n, 0, 0, 0, deparse.level = 0)

  ## If we never return a likelihood, the we never shuffle trajectories...
  p <- particle_filter$new(dat$data, dat$model, function(...) NULL)
  ll <- p$run(y0, n_particles, TRUE)
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
    if (observed$incidence > 15) {
      ret[] <- -Inf
    }
    ret
  }

  p <- particle_filter$new(dat$data, dat$model, compare)
  res <- p$run(dat$y0, 42, TRUE)
  expect_equal(res, -Inf)

  i <- (which(dat$data$incidence > 15)[[1]] + 2):101
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


test_that("predict", {
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  res <- p$run(dat$y0, 42, TRUE)
  t <- 0:10

  set.seed(1)
  res1 <- p$predict(t)
  set.seed(1)
  res2 <- p$predict(t, TRUE)

  expect_equal(dim(res1), c(5, 42, 10))
  expect_equal(dim(res2), c(5, 42, 111))

  ## history is prepended
  expect_equal(res2[, , 1:101], p$history)

  ## state is last in that
  expect_equal(res2[, , 101], p$state)

  ## appended predictions match the raw preductions
  expect_equal(res2[, , 102:111], res1)
})


test_that("can't predict until model has been run", {
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  expect_error(p$predict(0:10), "Particle filter has not been run")
})


test_that("can't append predictions without history", {
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  res <- p$run(dat$y0, 42, FALSE)
  expect_error(p$predict(0:10, TRUE), "Can't append without history")
})


test_that("prediction time must start at zero", {
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  res <- p$run(dat$y0, 42, FALSE)
  expect_error(p$predict(1:10), "Expected first 't' element to be zero")
})


test_that("Validate steps", {
  steps <- cbind(0:10 * 10, 1:11 * 10)
  expect_identical(particle_steps(steps, NULL), steps)
  expect_identical(particle_steps(steps, 0), steps)
  res <- particle_steps(steps, 5)
  expect_identical(res[-1], steps[-1])
  expect_identical(res[[1]], 5)

  expect_error(
    particle_steps(steps, -5),
    "'step_start' must be >= 0 (the first value of data$step_start)",
    fixed = TRUE)
  expect_error(
    particle_steps(steps, 10),
    "'step_start' must be < 10 (the first value of data$step_end)",
    fixed = TRUE)
  expect_error(
    particle_steps(steps, 20),
    "'step_start' must be < 10 (the first value of data$step_end)",
    fixed = TRUE)
})
