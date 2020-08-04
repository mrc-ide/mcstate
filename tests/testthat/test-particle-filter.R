context("particle_filter")

test_that("run particle filter on sir model", {
  dat <- example_sir()
  n_particles <- 42
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  res <- p$run()
  expect_is(res, "numeric")

  expect_is(p$state, "matrix")
  expect_equal(dim(p$state), c(3, n_particles))
  expect_equal(length(p$unique_particles), nrow(dat$data) + 1)
  expect_true(all(p$unique_particles <= n_particles & p$unique_particles >= 1))
  expect_null(p$history)
})


test_that("continuing a particle filter continues the RNG", {
  dat <- example_sir()
  n_particles <- 42
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  set.seed(1) # affects sample() used for filtering
  res <- p$run()
  expect_is(res, "numeric")

  set.seed(1)
  res2 <- p$run()
  expect_true(res2 != res)
})


test_that("run particle filter without index", {
  dat <- example_sir()
  n_particles <- 42
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)

  compare2 <- function(state, prev_state, ...) {
    dat$compare(state[4, , drop = FALSE], prev_state[4, , drop = FALSE], ...)
  }

  p2 <- particle_filter$new(dat$data, dat$model, n_particles, compare2)

  set.seed(1)
  ll1 <- p1$run()
  set.seed(1)
  ll2 <- p2$run()
  expect_identical(ll1, ll2)

  expect_equal(dim(p1$state), c(3, n_particles))
  expect_equal(dim(p2$state), c(4, n_particles))
})


test_that("particle filter likelihood is worse with worse parameters", {
  dat <- example_sir()
  n_particles <- 100
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  ll1 <- p$run()
  ll2 <- p$run(pars = list(gamma = 1, beta = 1))
  expect_true(ll1 > ll2)
})


test_that("stop simulation when likelihood is impossible", {
  dat <- example_sir()
  n_particles <- 42
  steps <- nrow(dat$data) + 1

  compare <- function(state, output, observed, pars) {
    ret <- dat$compare(state, output, observed, pars)
    if (observed$incidence > 15) {
      ret[] <- -Inf
    }
    ret
  }

  p <- particle_filter$new(dat$data, dat$model, n_particles, compare,
                           index = dat$index)
  res <- p$run(save_history = TRUE)
  expect_equal(res, -Inf)

  i <- (which(dat$data$incidence > 15)[[1]] + 2):steps
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
  n_particles <- 42
  steps <- nrow(dat$data) + 1

  set.seed(1)
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  run1 <- p1$run(save_history = TRUE)
  set.seed(1)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  run2 <- p2$run(save_history = TRUE)

  t <- 0:10
  res1 <- p1$predict(t)
  res2 <- p2$predict(t, TRUE)

  expect_equal(dim(res1), c(3, n_particles, length(t) - 1))
  expect_equal(dim(res2), c(3, n_particles, steps + length(t) - 1))

  ## history is prepended
  expect_equal(res2[, , 1:steps], p2$history)

  ## appended predictions match the raw predictions
  expect_equal(res2[, , (steps + 1):(steps + length(t) - 1)], res1)
})


test_that("predict particle filter without index", {
  dat <- example_sir()
  n_particles <- 42
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)

  compare2 <- function(state, prev_state, ...) {
    dat$compare(state[4, , drop = FALSE], prev_state[4, , drop = FALSE], ...)
  }

  p2 <- particle_filter$new(dat$data, dat$model, n_particles, compare2)

  set.seed(1)
  ll1 <- p1$run(save_history = TRUE)
  set.seed(1)
  ll2 <- p2$run(save_history = TRUE)
  expect_identical(ll1, ll2)

  t <- 0:10
  y1 <- p1$predict(t)
  y2 <- p2$predict(t)

  expect_equal(y1, y2[1:3, , ])
  expect_equal(y2[4, , ], 1000 - y1[1, , ])
})


test_that("can't predict until model has been run", {
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, 100, dat$compare,
                           index = dat$index)
  expect_error(p$predict(0:10), "Particle filter has not been run")
})


test_that("can't append predictions without history", {
  dat <- example_sir()
  n_particles <- 42
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  res <- p$run(save_history = FALSE)
  expect_error(p$predict(0:10, TRUE), "Can't append without history")
})


test_that("prediction time must start at zero", {
  dat <- example_sir()
  n_particles <- 42
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  res <- p$run(save_history = FALSE)
  expect_error(p$predict(1:10), "Expected first 't' element to be zero")
})


test_that("can't run predictions twice", {
  dat <- example_sir()
  n_particles <- 42
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  res <- p$run(save_history = FALSE)
  res <- p$predict(0:10)
  expect_error(p$predict(0:10), "Can't yet run predict multiple times")
})


test_that("Validate steps", {
  steps <- cbind(0:10 * 10, 1:11 * 10)
  expect_identical(particle_steps(steps, NULL), steps)
  expect_identical(particle_steps(steps, 0), steps)
  res <- particle_steps(steps, 5)
  expect_identical(res[-1], steps[-1])
  expect_identical(res[[1]], 5)

  res[1] <- 10L
  expect_identical(particle_steps(steps, 10), res)

  expect_error(
    particle_steps(steps, -5),
    "'step_start' must be >= 0 (the first value of data$step_start)",
    fixed = TRUE)
  expect_error(
    particle_steps(steps, 11),
    "'step_start' must be <= 10 (the first value of data$step_end)",
    fixed = TRUE)
  expect_error(
    particle_steps(steps, 20),
    "'step_start' must be <= 10 (the first value of data$step_end)",
    fixed = TRUE)
})


test_that("Control the comparison function", {
  dat <- example_sir()
  n_particles <- 42
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)

  pars <- list(compare = list(exp_noise = 1))
  ll1 <- p$run(pars = pars)

  pars$compare$exp_noise <- 0.01
  ll2 <- p$run(pars = pars)
  expect_true(ll2 < ll1)
})


test_that("Control the starting point of the simulation", {
  initial <- function(info, n_particles, pars) {
    list(step = pars$initial)
  }

  dat <- example_sir()
  n_particles <- 42
  data <- dat$data
  offset <- 400
  data[c("step_start", "step_end")] <-
    data[c("step_start", "step_end")] + offset
  data$step_start[[1]] <- 0

  ## The usual version:
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  set.seed(1)
  ll1 <- p1$run()

  ## Tuning the start date
  p2 <- particle_filter$new(data, dat$model, n_particles, dat$compare,
                            index = dat$index, initial = initial)
  set.seed(1)
  ll2 <- p2$run(list(initial = as.integer(offset)))
  expect_identical(ll1, ll2)

  ## Running from the beginning is much worse:
  set.seed(1)
  ll3 <- p2$run(list(initial = 0L))
  expect_true(ll3 < ll1)
})


test_that("control filter", {
  dat <- example_sir()
  expect_error(
    particle_filter$new(dat$data, dat$model, 0, dat$compare),
    "'n_particles' must be at least 1")
  expect_error(
    particle_filter$new(dat$data, dat$model, 1, dat$compare, seed = -1),
    "'seed' must be at least 1")
  expect_error(
    particle_filter$new(dat$data, dat$model, 1, dat$compare, n_threads = -1),
    "'n_threads' must be at least 1")
})


test_that("run particle filter on sir model", {
  dat <- example_sir()
  expect_error(
    particle_filter$new(dat$data, NULL, dat$compare, index = dat$index),
    "'model' must be a dust_generator")
})


test_that("scale log weights", {
  expect_equal(scale_log_weights(c(-Inf, -Inf)),
               list(weights = c(NaN, NaN), average = -Inf))
  expect_equal(scale_log_weights(c(-Inf, 1)),
               list(weights = c(0, 1), average = log(exp(1) / 2)))
  expect_equal(scale_log_weights(c(-Inf, 1, 1)),
               list(weights = c(0, 1, 1), average = log(exp(1) * 2 / 3)))
  expect_equal(scale_log_weights(c(NaN, NaN)),
               list(weights = c(NaN, NaN), average = -Inf))
  expect_equal(scale_log_weights(c(NaN, NaN)),
               list(weights = c(NaN, NaN), average = -Inf))
  expect_equal(scale_log_weights(c(NaN, 1)),
               list(weights = c(0, 1), average = log(exp(1) / 2)))
})


test_that("index must be sensible", {
  dat <- example_sir()
  expect_error(
    particle_filter$new(dat$data, dat$model, dat$compare,
                        index = c(1, 3, 5)),
    "'index' must be function if not NULL")
})


test_that("initial must be sensible", {
  dat <- example_sir()
  expect_error(
    particle_filter$new(dat$data, dat$model, dat$compare,
                        initial = c(1, 3, 5)),
    "'initial' must be function if not NULL")
})


test_that("we do not reorder particles when compare is NULL", {
  dat <- example_sir()
  n_particles <- 42
  p <- particle_filter$new(dat$data, dat$model, n_particles,
                           function(...) NULL, index = dat$index)
  res <- p$run()
  expect_equal(res, 0)
  expect_equal(p$unique_particles, rep(n_particles, 101))
})


test_that("initialise with simple state", {
  dat <- example_sir()
  n_particles <- 100

  initial <- function(info, n_particles, pars) {
    c(1000, pars$I0, 0, 0)
  }
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, initial = initial)

  ll1 <- p$run(list(I0 = 200), save_history = TRUE)
  expect_equal(p$history[, , 1],
               matrix(c(1000, 200, 0), 3, n_particles))
  ll2 <- p$run(list(I0 = 1), save_history = TRUE)
  expect_equal(p$history[, , 1],
               matrix(c(1000, 1, 0), 3, n_particles))
  ll3 <- p$run(list(I0 = 10), save_history = TRUE)
  expect_equal(p$history[, , 1],
               matrix(c(1000, 10, 0), 3, n_particles))

  expect_true(ll1 < ll3)
  expect_true(ll2 < ll3)
})


test_that("initialise with complex state", {
  dat <- example_sir()
  n_particles <- 100

  initial <- function(info, n_particles, pars) {
    y <- matrix(0, 4, n_particles)
    set.seed(1) # so that we can check below
    i0 <- rpois(n_particles, pars$I0)
    y[1, ] <- 1100 - i0
    y[2, ] <- i0
    y
  }

  ## Set the incidence to NA so that no shuffling occurs
  dat$data$incidence <- NA
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, initial = initial)

  pars <- list(I0 = 10)

  set.seed(1)
  ll <- p$run(pars, save_history = TRUE)
  expect_equal(p$history[, , 1],
               initial(NULL, n_particles, pars)[1:3, ])
})


test_that("Control the starting point of the simulation", {
  dat <- example_sir()
  n_particles <- 42

  data <- dat$data
  ## Drop extra columns that are confusing in this context
  data <- data[setdiff(names(data), c("day_start", "day_end"))]
  ## Offset the data by a considerable margin:
  offset <- 400
  data[c("step_start", "step_end")] <-
    data[c("step_start", "step_end")] + offset
  ## A burnin step:
  intro <- data[1, ]
  intro$incidence <- NA
  intro$step_start <- 0
  intro$step_end <- data$step_start[[1]]
  ## Our combination data:
  data <- rbind(intro, data)

  ## Now we can start work with an "interesting" starting point.

  ## The usual version, with normal data:
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  set.seed(1)
  ll1 <- p1$run()

  ## Then tune the start date to get the same effect:
  initial <- function(info, n_particles, pars) {
    list(state = c(1000, pars$I0, 0, 0),
         step = pars$step)
  }

  ## Tuning the start date
  p2 <- particle_filter$new(data, dat$model, n_particles, dat$compare,
                            index = dat$index, initial = initial)
  set.seed(1)
  ll2 <- p2$run(list(I0 = 10, step = offset),
                save_history = TRUE)

  expect_identical(ll1, ll2)

  ## Running from the beginning is much worse:
  set.seed(1)
  ll3 <- p2$run(list(I0 = 10, step = 0),
                save_history = TRUE)
  expect_true(ll3 < ll1)
})


test_that("Variable initial starting point of the simulation", {
  dat <- example_sir()
  n_particles <- 42

  data <- dat$data
  ## Drop extra columns that are confusing in this context
  data <- data[setdiff(names(data), c("day_start", "day_end"))]
  ## Offset the data by a considerable margin:
  offset <- 400
  data[c("step_start", "step_end")] <-
    data[c("step_start", "step_end")] + offset
  ## A burnin step:
  intro <- data[1, ]
  intro$incidence <- NA
  intro$step_start <- 0
  intro$step_end <- data$step_start[[1]]
  ## Our combination data:
  data <- rbind(intro, data)

  ## Then tune the start date to get the same effect:
  initial <- function(info, n_particles, pars) {
    set.seed(1)
    i0 <- rpois(n_particles, pars$I0)
    step <- pars$step_offset - rpois(n_particles, pars$step_mean)
    list(state = rbind(1000, i0, 0, 0, deparse.level = 0),
         step = step)
  }
  compare <- function(...) {
    NULL
  }

  p <- particle_filter$new(data, dat$model, n_particles, compare,
                           index = dat$index, initial = initial)
  pars <- list(I0 = 10, step_mean = 20, step_offset = 400)
  set.seed(1)
  ll <- p$run(pars, save_history = TRUE)
  expect_equal(ll, 0)

  mod <- p$model$new(list(), step = 0, n_particles = n_particles)
  tmp <- initial(NULL, n_particles, pars)
  mod$set_state(tmp$state, tmp$step)
  expect_equal(p$history[, , 1], mod$state()[1:3, ])
  expect_equal(mod$step(), max(tmp$step))
})
