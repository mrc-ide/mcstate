context("particle_filter")

test_that("run particle filter on sir model", {
  dat <- example_sir()
  n_particles <- 42
  set.seed(1)
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  res <- p$run()
  expect_is(res, "numeric")

  state <- p$state()
  expect_is(state, "matrix")
  expect_equal(dim(state), c(5, n_particles))

  expect_error(
    p$history(),
    "Can't get history as model was run with save_history = FALSE")
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

  compare2 <- function(state, ...) {
    dat$compare(state[5, , drop = FALSE], ...)
  }

  p2 <- particle_filter$new(dat$data, dat$model, n_particles, compare2)

  set.seed(1)
  ll1 <- p1$run(save_history = TRUE)
  set.seed(1)
  ll2 <- p2$run(save_history = TRUE)
  expect_identical(ll1, ll2)

  expect_equal(dim(p1$history()), c(3, n_particles, 101))
  expect_equal(dim(p2$history()), c(5, n_particles, 101))
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

  compare <- function(state, observed, pars) {
    ret <- dat$compare(state, observed, pars)
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
  history <- p$history()
  expect_false(any(is.na(history[, , !i])))
  expect_true(all(is.na(history[, , i])))
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
  dat <- example_sir()

  ## The usual version:
  n_particles <- 42
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  set.seed(1)
  ll1 <- p1$run()

  ## Tuning the start date
  data_raw <- dat$data_raw
  data_raw$day <- data_raw$day + 100
  data <- particle_filter_data(data_raw, "day", 4, 100)

  p2 <- particle_filter$new(data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  set.seed(1)
  ll2 <- p2$run()
  expect_identical(ll1, ll2)

  ## Running from the beginning is much worse:
  set.seed(1)
  data <- particle_filter_data(data_raw, "day", 4, 1)
  p3 <- particle_filter$new(data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  ll3 <- p3$run()
  expect_true(ll3 < ll1)
})


test_that("Cannot use previous initial condition approach", {
  initial <- function(info, n_particles, pars) {
    list(step = 2)
  }

  dat <- example_sir()
  n_particles <- 42
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, initial = initial)
  expect_error(p$run(), "Setting 'step' from initial no longer supported")
})


test_that("control filter", {
  dat <- example_sir()
  expect_error(
    particle_filter$new(dat$data, dat$model, 0, dat$compare),
    "'n_particles' must be at least 1")
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
})


test_that("initialise with simple state", {
  dat <- example_sir()
  n_particles <- 100

  initial <- function(info, n_particles, pars) {
    c(1000, pars$I0, 0, 0, 0)
  }
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, initial = initial)

  ll1 <- p$run(list(I0 = 200), save_history = TRUE)
  expect_equal(p$history()[, , 1],
               matrix(c(1000, 200, 0), 3, n_particles))
  ll2 <- p$run(list(I0 = 1), save_history = TRUE)
  expect_equal(p$history()[, , 1],
               matrix(c(1000, 1, 0), 3, n_particles))
  ll3 <- p$run(list(I0 = 10), save_history = TRUE)
  expect_equal(p$history()[, , 1],
               matrix(c(1000, 10, 0), 3, n_particles))

  expect_true(ll1 < ll3)
  expect_true(ll2 < ll3)
})


test_that("initialise with complex state", {
  dat <- example_sir()
  n_particles <- 100

  initial <- function(info, n_particles, pars) {
    y <- matrix(0, 5, n_particles)
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
  expect_equal(p$history()[, , 1],
               initial(NULL, n_particles, pars)[1:3, ])
})


test_that("can save history", {
  dat <- example_sir()
  n_particles <- 42
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  p$run(save_history = TRUE)
  res <- p$history()

  ## If we have correctly sampled trajectories, then we'll have
  ## monotonic S and R within a particle:
  expect_true(all(diff(t(res[1, , ])) <= 0))
  expect_true(all(diff(t(res[3, , ])) >= 0))

  ## Can get just a few histories
  expect_equal(
    drop(p$history(1)),
    res[, 1, ])
  expect_equal(
    drop(p$history(10)),
    res[, 10, ])
  expect_equal(
    drop(p$history(10:20)),
    res[, 10:20, ])
})


test_that("can't get state or history until model is run", {
  dat <- example_sir()
  n_particles <- 42
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  expect_error(
    p$state(),
    "Model has not yet been run")
  expect_error(
    p$history(),
    "Model has not yet been run")
  expect_error(
    p$restart_state(),
    "Model has not yet been run")
})


test_that("can filter state on extraction", {
  dat <- example_sir()
  n_particles <- 42
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  res <- p$run()
  expect_is(res, "numeric")

  state <- p$state()
  expect_equal(p$state(1), state[1, , drop = FALSE])
  state <- p$state()
  expect_equal(p$state(2:3), state[2:3, , drop = FALSE])
})


test_that("can return inputs", {
  dat <- example_sir()
  n_particles <- 42
  initial <- function(...) NULL
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, initial = initial, seed = 100)
  inputs <- p$inputs()
  expect_setequal(names(inputs), names(formals(p$initialize)))

  expect_equal(inputs$data, dat$data)
  expect_equal(inputs$model, dat$model)
  expect_equal(inputs$n_particles, n_particles)
  expect_equal(inputs$index, dat$index)
  expect_equal(inputs$compare, dat$compare)
  expect_null(inputs$gpu_config)
  expect_equal(inputs$initial, initial)
  expect_equal(inputs$seed, 100)

  res <- p$run()

  inputs2 <- p$inputs()
  expect_type(inputs2$seed, "raw")

  expect_identical(inputs2[names(inputs2) != "seed"],
                   inputs[names(inputs) != "seed"])
})


test_that("return names on history, if present", {
  dat <- example_sir()
  n_particles <- 42
  index <- function(info) {
    list(run = 4L, state = c(S = 1L, I = 2L, R = 3L))
  }
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = index, seed = 100)
  p$run(save_history = TRUE)
  res <- p$history()
  expect_equal(rownames(res), c("S", "I", "R"))
})


test_that("no names on history, if absent", {
  dat <- example_sir()
  n_particles <- 42
  index <- function(info) {
    list(run = 4L, state = 1:3)
  }
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = index, seed = 100)
  p$run(save_history = TRUE)
  res <- p$history()
  expect_null(rownames(res))
})


test_that("can change the number of threads (null model)", {
  skip_on_cran()
  dat <- example_sir()
  n_particles <- 42
  index <- function(info) {
    list(run = 4L, state = 1:3)
  }
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = index, seed = 100, n_threads = 1L)
  expect_equal(p$set_n_threads(2), 1)
  expect_equal(r6_private(p)$n_threads, 2)
  p$run()
  expect_equal(r6_private(p)$last_model$n_threads(), 2)

  expect_equal(p$set_n_threads(1), 2)
  expect_equal(r6_private(p)$n_threads, 1)
  expect_equal(r6_private(p)$last_model$n_threads(), 1)
  p$run()
  expect_equal(r6_private(p)$last_model$n_threads(), 1)
})


test_that("Can extract state from the model", {
  dat <- example_sir()
  n_particles <- 42
  seed <- 100
  index <- function(info) {
    list(run = 5L, state = 1:3)
  }
  set.seed(1)
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = index, seed = seed)
  end <- c(20, 40, 60)
  res <- p$run(save_restart = end)

  s <- p$restart_state()
  expect_equal(dim(s), c(5, n_particles, length(end)))

  f <- function(n) {
    set.seed(1)
    d <- dat$data[dat$data$day_end <= n, ]
    p <- particle_filter$new(d, dat$model, n_particles, dat$compare,
                             index = index, seed = seed)
    p$run()
    p$state()
  }

  cmp <- lapply(end, f)
  expect_equal(s[, , 1], cmp[[1]])
  expect_equal(s[, , 2], cmp[[2]])
  expect_equal(s[, , 3], cmp[[3]])

  expect_equal(p$restart_state(1), s[, 1, , drop = FALSE])
  expect_equal(p$restart_state(c(10, 3, 3, 6)), s[, c(10, 3, 3, 6), ])
})


test_that("can extract just one restart state", {
  dat <- example_sir()
  n_particles <- 42
  seed <- 100
  index <- function(info) {
    list(run = 5L, state = 1:3)
  }
  set.seed(1)
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = index, seed = seed)
  end <- 40
  res <- p$run(save_restart = end)
  s <- p$restart_state()
  s1 <- p$restart_state(1)
  expect_equal(dim(s), c(5, n_particles, 1))
  expect_equal(dim(s1), c(5, 1, 1))

  set.seed(1)
  d <- dat$data[dat$data$day_end <= end, ]
  p <- particle_filter$new(d, dat$model, n_particles, dat$compare,
                           index = index, seed = seed)
  p$run()
  cmp <- p$state()
  expect_equal(drop(s), cmp)
  expect_equal(drop(s1), cmp[, 1])
})


test_that("Can't get restart state without saving it", {
  dat <- example_sir()
  index <- function(info) {
    list(run = 4L, state = 1:3)
  }
  set.seed(1)
  p <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                           index = index, seed = 100)
  res <- p$run()
  expect_error(p$restart_state(),
               "Can't get history as model was run with save_restart = NULL")
})


test_that("Dates must exist to save restart state", {
  dat <- example_sir()
  index <- function(info) {
    list(run = 4L, state = 1:3)
  }
  set.seed(1)
  p <- particle_filter$new(dat$data[1:20, ], dat$model, 42, dat$compare,
                           index = index, seed = 100)
  expect_error(
    p$run(save_restart = 30),
    "'save_restart' contains times not in 'day': 30")
  expect_error(
    p$run(save_restart = c(30, 50)),
    "'save_restart' contains times not in 'day': 30, 50")
  expect_error(
    p$run(save_restart = c(10, 30, 50)),
    "'save_restart' contains times not in 'day': 30, 50")
})


test_that("use compiled compare function", {
  dat <- example_sir()
  n_particles <- 100
  set.seed(1)

  model <- dust::dust_example("sir")
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, model, n_particles, NULL,
                            index = dat$index)

  ## Proving these are the same is tricky to do in a sensible amount
  ## of time but if we run 1000 replicates with 400 particles it's
  ## easy to see that these are the same (this just takes the best
  ## part of a minute and replicates the unit tests available
  ## elsewhere)
  y1 <- replicate(50, p1$run())
  y2 <- replicate(50, p2$run())
  expect_equal(mean(y1), mean(y2), tolerance = 0.01)
})


test_that("Can get history with compiled particle filter", {
  dat <- example_sir()
  n_particles <- 100
  set.seed(1)

  model <- dust::dust_example("sir")
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, model, n_particles, NULL,
                            index = dat$index)

  p1$run(save_history = TRUE)
  p2$run(save_history = TRUE)

  expect_equal(dim(p1$history()), dim(p2$history()))
  expect_true(all(diff(t(p2$history()[3, , ])) >= 0))
  expect_equal(dim(p1$history(1L)), dim(p2$history(1L)))
  expect_equal(dim(p1$history(1:5)), dim(p2$history(1:5)))
})


test_that("Can save restart with compiled particle filter", {
  dat <- example_sir()
  n_particles <- 100
  set.seed(1)

  model <- dust::dust_example("sir")
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, model, n_particles, NULL,
                            index = dat$index)

  at <- c(10, 20, 40, 80)
  p1$run(save_restart = at)
  p2$run(save_restart = at)

  expect_equal(dim(p1$restart_state()), dim(p2$restart_state()))
  expect_true(all(diff(t(p2$restart_state()[3, , ])) >= 0))
  expect_equal(dim(p1$restart_state(1L)), dim(p2$restart_state(1L)))
  expect_equal(dim(p1$restart_state(1:5)), dim(p2$restart_state(1:5)))
})


test_that("prevent using compiled compare where model does not support it", {
  dat <- example_sir()
  n_particles <- 100
  model <- dust::dust_example("walk")
  expect_error(
    particle_filter$new(dat$data, model, n_particles, NULL),
    "Your model does not have a built-in 'compare' function")
})


test_that("can't get partial likelihood with compiled compare", {
  dat <- example_sir()
  n_particles <- 100
  set.seed(1)

  model <- dust::dust_example("sir")
  p <- particle_filter$new(dat$data, model, n_particles, NULL,
                           index = dat$index)
  obj <- p$run_begin()
  expect_error(
    obj$step(10, TRUE),
    "'partial' not supported with compiled compare")
})


test_that("incrementally run a particle filter", {
  dat <- example_sir()
  n_particles <- 42

  set.seed(1)
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  cmp <- p1$run()

  set.seed(1)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)

  n <- nrow(dat$data)
  ans1 <- numeric(n)
  ans2 <- numeric(n)
  obj <- p2$run_begin()
  expect_s3_class(obj, "particle_filter_state")
  for (i in seq_len(n)) {
    ans1[[i]] <- obj$step(i)
    ans2[[i]] <- obj$log_likelihood_step
  }

  expect_identical(obj$log_likelihood, cmp)
  expect_identical(ans1[[length(ans1)]], cmp)
  expect_true(all(diff(ans1) < 0))
  expect_equal(c(ans1[[1]], diff(ans1)), ans2)
})


test_that("incrementally run a compiled particle filter", {
  dat <- example_sir()
  n_particles <- 42

  set.seed(1)
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, NULL,
                           index = dat$index, seed = 1L)
  cmp <- p1$run()

  set.seed(1)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, NULL,
                            index = dat$index, seed = 1L)

  n <- nrow(dat$data)
  ans <- numeric(n)
  obj <- p2$run_begin()
  expect_s3_class(obj, "particle_filter_state")
  for (i in seq_len(n)) {
    ans[[i]] <- obj$step(i)
  }

  expect_identical(obj$log_likelihood, cmp)
  expect_identical(ans[[length(ans)]], cmp)
  expect_true(all(diff(ans) < 0))
})


test_that("Can't step a particle filter object past its end", {
  dat <- example_sir()
  n_particles <- 42

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  obj <- p$run_begin()
  n <- nrow(dat$data)
  obj$run()
  expect_error(obj$run(),
               "Particle filter has reached the end of the data")
  expect_error(obj$step(n),
               "Particle filter has reached the end of the data")
  expect_error(obj$step(n + 1),
               "Particle filter has reached the end of the data")
})


test_that("Can't rerun a step", {
  dat <- example_sir()
  n_particles <- 42

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  obj <- p$run_begin()
  n <- nrow(dat$data)
  obj$step(10)
  expect_error(
    obj$step(10),
    "Particle filter has already run step index 10 (to model step 40)",
    fixed = TRUE)
  expect_error(
    obj$step(5),
    "Particle filter has already run step index 5 (to model step 20)",
    fixed = TRUE)
})


test_that("Can't run past the end of the data", {
  dat <- example_sir()
  n_particles <- 42

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  obj <- p$run_begin()
  n <- nrow(dat$data)
  expect_error(
    obj$step(n + 1),
    "step_index 101 is beyond the length of the data (max 100)",
    fixed = TRUE)
})


test_that("Can fork a particle_filter_state object", {
  dat <- example_sir()
  n_particles <- 42

  set.seed(1)
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  ans <- numeric(nrow(dat$data))
  obj <- p1$run_begin(save_history = TRUE)

  for (i in seq_len(10)) {
    ans[[i]] <- obj$step(i)
  }

  tmp <- obj$model$rng_state()

  set.seed(1)
  res <- obj$fork_smc2(list(beta = 0.15))
  expect_identical(res$model$rng_state(), obj$model$rng_state())
  expect_false(identical(obj$model$rng_state(), tmp))

  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = tmp)
  set.seed(1)
  cmp <- p2$run_begin(list(beta = 0.15), save_history = TRUE)
  cmp$step(10)

  expect_identical(res$log_likelihood, cmp$log_likelihood)
  expect_identical(res$history, cmp$history)
})

test_that("run particle filter on shared sir model", {
  dat <- example_sir_shared()
  n_particles <- 42

  pars <- list(list(beta = 0.2, gamma = 0.1),
                               list(beta = 0.3, gamma = 0.1))

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  set.seed(1)
  res <- p$run(pars)
  expect_is(res, "numeric")
  expect_equal(length(res), 2)

  state <- p$state()
  expect_is(state, "array")
  expect_equal(dim(state), c(5, n_particles, 2))

  expect_error(
    p$history(),
    "Can't get history as model was run with save_history = FALSE")

  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  set.seed(1)
  expect_identical(res, p2$run(pars))
})

test_that("can save history - nested", {
  dat <- example_sir_shared()
  n_particles <- 42
  set.seed(1)

  pars <- list(list(beta = 0.2, gamma = 0.1),
                               list(beta = 0.3, gamma = 0.1))

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  p$run(pars, save_history = TRUE)
  res <- p$history()

  expect_equal(dim(res), c(3, 42, 2, 101))
  expect_equal(dim(p$history(1)), c(3, 1, 2, 101))
  expect_error(p$history(matrix(1, ncol = 3)), "2 columns")

  ## If we have correctly sampled trajectories, then we'll have
  ## monotonic S and R within a particle:
  expect_true(all(diff(t(res[1, , 1, ])) <= 0))
  expect_true(all(diff(t(res[1, , 2, ])) <= 0))
  expect_true(all(diff(t(res[3, , 1, ])) >= 0))
  expect_true(all(diff(t(res[3, , 2, ])) >= 0))

  ## Can get just a few histories
  expect_equal(
    drop(p$history(1)),
    res[, 1, , ])
  expect_equal(
    drop(p$history(10)),
    res[, 10, , ])
  expect_equal(
    drop(p$history(10:20)),
    res[, 10:20, , ])
})

test_that("Can extract state from the model - nested", {
  dat <- example_sir_shared()
  n_particles <- 42
  set.seed(1)
  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.1))
  seed <- 100
  index <- function(info) {
    list(run = 5L, state = 1:3)
  }

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = seed)

  end <- c(20, 40, 60)
  res <- p$run(pars, save_restart = end)

  s <- p$restart_state()
  expect_equal(dim(s), c(5, n_particles, 2, length(end)))

  f <- function(n) {
    set.seed(1)
    d <- dat$data[dat$data$day_end <= n, ]
    p <- particle_filter$new(d, dat$model, n_particles, dat$compare,
                             index = index, seed = seed)
    p$run(pars)
    p$state()
  }

  cmp <- lapply(end, f)
  expect_equal(s[, , , 1], cmp[[1]])
  expect_equal(s[, , , 2], cmp[[2]])
  expect_equal(s[, , , 3], cmp[[3]])

  expect_equal(p$restart_state(1), s[, 1, , , drop = FALSE])
  expect_equal(p$restart_state(c(10, 3, 3, 6)), s[, c(10, 3, 3, 6), , ])
})


test_that("use compiled compare function - nested", {
  dat <- example_sir_shared()
  n_particles <- 42
  set.seed(1)

  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.1))

  model <- dust::dust_example("sir")
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, model, n_particles, NULL,
                            index = dat$index)

  y1 <- replicate(50, p1$run(pars))
  y2 <- replicate(50, p2$run(pars))
  expect_equal(mean(y1), mean(y2), tolerance = 0.01)
})


test_that("can get history with compiled particle filter on nested model", {
  dat <- example_sir_shared()
  n_particles <- 42
  set.seed(1)

  pars <- list(list(beta = 0.2, gamma = 0.1),
                               list(beta = 0.3, gamma = 0.1))

  model <- dust::dust_example("sir")
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  p2 <- particle_filter$new(dat$data, model, n_particles, NULL,
                            index = dat$index)

  p1$run(pars, save_history = TRUE)
  p2$run(pars, save_history = TRUE)

  expect_equal(dim(p1$history()), dim(p2$history()))
  expect_true(all(diff(t(p2$history()[3, , 1, ])) >= 0))
  expect_true(all(diff(t(p2$history()[3, , 2, ])) >= 0))
  expect_equal(dim(p1$history(1L)), dim(p2$history(1L)))
  expect_equal(dim(p1$history(1:5)), dim(p2$history(1:5)))

  idx <- cbind(1:4, 2:5)
  h2 <- p2$history(idx)
  expect_equal(dim(h2), dim(p1$history(idx)))
  expect_equal(h2[, , 1, ], p2$history()[, 1:4, 1, ])
  expect_equal(h2[, , 2, ], p2$history()[, 2:5, 2, ])

  expect_error(p2$history(idx[, c(1, 1, 2)]),
               "'index_particle' should have 2 columns")
})


test_that("particle filter state nested - errors", {
  dat <- example_sir_shared()
  n_particles <- 42
  set.seed(1)

  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.1))

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)

  expect_error(p$run(pars[1]),
               "'pars' must have length 2",
               fixed = TRUE)
})

test_that("can't change initial step via initial in nested filter", {
  initial <- function(info, n_particles, pars) {
    list(step = 0)
  }

  dat <- example_sir_shared()
  n_particles <- 42

  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.1))

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, initial = initial)
  expect_error(p$run(pars),
               "Setting 'step' from initial no longer supported")
})

test_that("Can fork a particle_filter_state_nested object", {
  dat <- example_sir_shared()
  n_particles <- 42

  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.1))

  set.seed(1)
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)

  ans <- matrix(nrow = nrow(dat$data), ncol = 2)
  obj <- p1$run_begin(pars, save_history = TRUE)

  for (i in seq_len(10)) {
    ans[i, ] <- obj$step(i)
  }

  tmp <- obj$model$rng_state()

  set.seed(1)
  res <- obj$fork_smc2(pars)
  expect_identical(res$model$rng_state(), obj$model$rng_state())
  expect_false(identical(obj$model$rng_state(), tmp))

  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = tmp)
  set.seed(1)
  cmp <- p2$run_begin(pars, save_history = TRUE)
  cmp$step(10)

  expect_identical(res$log_likelihood, cmp$log_likelihood)
  expect_identical(res$history, cmp$history)
})

test_that("particle filter state nested - error steps", {
  dat <- example_sir_shared()
  n_particles <- 42

  pars <- list(list(beta = 0.2, gamma = 0.1),
                               list(beta = 0.3, gamma = 0.1))

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  set.seed(1)
  obj <- p$run_begin(pars)
  n <- nrow(dat$data)
  obj$run()
  expect_error(obj$run(),
               "Particle filter has reached the end of the data")
  expect_error(obj$step(n),
               "Particle filter has reached the end of the data")
  expect_error(obj$step(n + 1),
               "Particle filter has reached the end of the data")

  obj <- p$run_begin(pars)
  obj$step(10)
  expect_error(
    obj$step(10),
    "Particle filter has already run step index 10 (to model step 40)",
    fixed = TRUE)
  expect_error(
    obj$step(5),
    "Particle filter has already run step index 5 (to model step 20)",
    fixed = TRUE)

  expect_error(
    obj$step(n + 1),
    "step_index 201 is beyond the length of the data (max 100)",
    fixed = TRUE)
})

test_that("stop simulation when likelihood is impossible", {
  dat <- example_sir_shared()
  n_particles <- 42

  pars <- list(list(beta = 0.2, gamma = 0.1),
                               list(beta = 0.3, gamma = 0.1))

  n_particles <- 42
  steps <- nrow(dat$data) / 2 + 1

  compare <- function(state, observed, pars) {
    ret <- dat$compare(state, observed, pars)
    if (observed$incidence > 15) {
      ret[] <- -Inf
    }
    ret
  }

  p <- particle_filter$new(dat$data, dat$model, n_particles, compare,
                           index = dat$index)
  res <- p$run(pars, save_history = TRUE)
  expect_true(-Inf %in% res)

  i <- (which(dat$data$incidence[1:100] > 15)[[1]] + 2):steps
  history <- p$history()
  expect_false(any(is.na(history[, , 1, !i])))
  expect_true(all(is.na(history[, , 1, i])))
})

test_that("compare NULL - nested", {
  dat <- example_sir_shared()
  n_particles <- 42

  pars <- list(list(beta = 0.2, gamma = 0.1),
                               list(beta = 0.3, gamma = 0.1))
  p <- particle_filter$new(dat$data, dat$model, n_particles,
                           function(...) NULL, index = dat$index)
  res <- p$run(pars, save_history = TRUE)
  expect_equal(res, c(0, 0))
  expect_equal(dim(p$history()), c(3, 42, 2, 101))
})

test_that("nested particle filter initial not list", {
  dat <- example_sir_shared()
  n_particles <- 42
  initial <- function(...) NULL
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, initial = initial, seed = 100)
  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.1))
  expect_is(p$run(pars), "numeric")
})


test_that("return names on nested history, if present", {
  dat <- example_sir_shared()
  n_particles <- 42
  index <- function(info) {
    list(run = 4L, state = c(S = 1L, I = 2L, R = 3L))
  }
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
    index = index, seed = 100
  )
  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.1))
  p$run(pars, save_history = TRUE)
  res <- p$history()
  expect_equal(rownames(res), c("S", "I", "R"))
})


test_that("error on different population indices", {
  dat <- example_sir_shared()
  n_particles <- 42
  index <- function(info) {
    list(run = info$pars$beta * 10)
  }
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
    index = index, seed = 100
  )
  pars <- list(
    list(beta = 0.2, gamma = 0.1),
    list(beta = 0.3, gamma = 0.1)
  )
  set.seed(1)
  expect_error(p$run(pars, save_history = TRUE),
               "index must be identical across populations")
})

test_that("initialise with complex state - nested", {
  dat <- example_sir_shared()
  n_particles <- 100

  initial <- function(info, n_particles, pars) {
    y <- matrix(0, 5, n_particles)
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

  pars <- list(list(I0 = 10, beta = 0.2, gamma = 0.1),
               list(I0 = 10, beta = 0.3, gamma = 0.1))

  set.seed(1)
  ll <- p$run(pars, save_history = TRUE)
  expect_equal(p$history()[, , 1, 1],
               initial(NULL, n_particles, pars[[1]])[1:3, ])
  expect_equal(p$history()[, , 2, 1],
               initial(NULL, n_particles, pars[[2]])[1:3, ])
})


test_that("Control the starting point of a nested simulation", {
  dat <- example_sir_shared()
  n_particles <- 42

  pars <- list(list(I0 = 10, beta = 0.2, gamma = 0.1),
               list(I0 = 10, beta = 0.3, gamma = 0.1))

  ## The usual version
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  set.seed(1)
  ll1 <- p1$run(pars)

  ## Tuning the start date
  data_raw <- dat$data_raw
  data_raw$day <- data_raw$day + 100
  data <- particle_filter_data(data_raw, time = "day", 4, 100,
                               population = "populations")

  p2 <- particle_filter$new(data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  set.seed(1)
  ll2 <- p2$run(pars)
  expect_identical(ll1, ll2)
})


test_that("nested error on unequal state", {
  dat <- example_sir_shared()
  n_particles <- 42
  pars <- list(list(I0 = NULL, beta = 0.2, gamma = 0.1),
               list(I0 = 10, beta = 0.3, gamma = 0.1))

  initial <- function(info, n_particles, pars) {
    c(1000, pars$I0, 0, 0, 0)
  }

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, initial = initial)

  expect_error(p$run(pars), "unequal state")
})


test_that("nested silent on initial w. state w/o step", {
  dat <- example_sir_shared()
  n_particles <- 42
  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.1))

  initial <- function(info, n_particles, pars) {
    c(1000, 0, 0, 0, 0)
  }

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, initial = initial)
  expect_is(p$run(pars), "numeric")
})


test_that("gpu_config requires GPU support", {
  dat <- example_sir()
  n_particles <- 100
  set.seed(1)
  expect_error(
    particle_filter$new(dat$data, dat$model, n_particles, NULL,
                        index = dat$index, gpu_config = 0),
    "gpu_config' provided, but 'model' does not have GPU support")
})


test_that("Can run a gpu model by passing gpu_config through", {
  dat <- example_sir()
  n_particles <- 100
  set.seed(1)

  model <- dust::dust_example("sirs")
  p_c <- particle_filter$new(dat$data, model, n_particles, NULL,
                             index = dat$index)
  p_g <- particle_filter$new(dat$data, model, n_particles, NULL,
                             index = dat$index, gpu_config = 0)
  expect_null(r6_private(p_c)$gpu_config)
  expect_equal(r6_private(p_g)$gpu_config, 0)

  filter_c <- p_c$run_begin()
  filter_g <- p_g$run_begin()

  expect_false(r6_private(filter_c)$gpu)
  expect_true(r6_private(filter_g)$gpu)
  expect_false(filter_c$model$uses_gpu(TRUE))
  expect_true(filter_g$model$uses_gpu(TRUE))
})


test_that("Can terminate a filter early", {
  dat <- example_sir()
  n_particles <- 42
  ## First run through, we'll get our cutoffs:
  set.seed(1)
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  ll1 <- replicate(10, p1$run())

  set.seed(1)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  min_ll <- mean(ll1)
  ll2 <- replicate(10, p2$run(min_log_likelihood = min_ll))
  expect_true(-Inf %in% ll2)
  expect_true(min(ll2[is.finite(ll2)]) >= min_ll)
})


test_that("nested particle filter requires unnamed parameters", {
  dat <- example_sir_shared()
  p <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                           index = dat$index)
  pars <- list(a = list(beta = 0.2, gamma = 0.1),
               b = list(beta = 0.3, gamma = 0.1))
  expect_error(
    p$run(pars),
    "Expected an unnamed list of parameters")
})


test_that("Can do early exit for nested filter", {
  dat <- example_sir_shared()
  n_particles <- 42
  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.1))

  ## Same setup as before
  set.seed(1)
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  ll1 <- replicate(10, p1$run(pars))

  set.seed(1)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  min_ll <- mean(colSums(ll1))
  ll2 <- replicate(10, p2$run(pars, min_log_likelihood = min_ll))
  expect_true(-Inf %in% ll2)
  expect_true(min(ll2[is.finite(ll2)]) >= min_ll)

  set.seed(1)
  p3 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  min_ll <- ll1[, which.min(abs(colSums(ll1) - mean(colSums(ll1))))]
  ll3 <- replicate(10, p3$run(pars, min_log_likelihood = min_ll))

  i <- apply(ll3 > -Inf, 2, any)
  expect_true(!all(i))
  expect_true(all(ll3[, !i] == -Inf))
  expect_true(all(apply(ll3[, i] >= min_ll, 2, any)))
  expect_false(all(apply(ll3[, i] >= min_ll, 2, all)))
})


test_that("Confirm nested filter is correct", {
  ## To show this works, we'll run the filters separately.  This will
  ## be a useful result when adding multistage parameters.

  ## This example is a bit fiddly to get exact equivalence - we need
  ## to manually step the filter to get uses of R's RNG to be correct.
  ## This would go away if we updated to use the same strategy as the
  ## the compiled filter and use a dust RNG here, stepping off the end
  ## of the last rng state.  That requires some additional work though
  ## to keep the state in sync.
  ##
  ## The other bit of fiddle is replacing the stochastic noise in
  ## compare with deterministic epsilon to avoid -Inf likelihoods in
  ## dpois
  dat <- example_sir_shared()
  n_particles <- 42

  ## The usual compare, but add a fixed amount of noise
  compare <- function(state, observed, pars = NULL) {
    if (is.na(observed$incidence)) {
      return(NULL)
    }
    incidence_modelled <- state[1, , drop = TRUE]
    incidence_observed <- observed$incidence
    lambda <- incidence_modelled + 1e-7
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
  }

  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.1))

  seed <- dust::dust_rng_pointer$new(1, n_particles * 2)$state()
  seed1 <- seed[seq_len(length(seed) / 2)]
  seed2 <- seed[-seq_len(length(seed) / 2)]

  data1 <- particle_filter_data(dat$data_raw[dat$data_raw$populations == "a", ],
                                time = "day", rate = 4)
  data2 <- particle_filter_data(dat$data_raw[dat$data_raw$populations == "b", ],
                                time = "day", rate = 4)

  set.seed(1)
  p1 <- particle_filter$new(data1, dat$model, n_particles, compare,
                            index = dat$index, seed = seed1)
  p2 <- particle_filter$new(data2, dat$model, n_particles, compare,
                            index = dat$index, seed = seed2)

  s1 <- p1$run_begin(pars[[1]], save_history = TRUE)
  s2 <- p2$run_begin(pars[[2]], save_history = TRUE)
  for (i in seq_len(nrow(data1))) {
    s1$step(i)
    s2$step(i)
  }

  set.seed(1)
  p3 <- particle_filter$new(dat$data, dat$model, n_particles, compare,
                            index = dat$index, seed = seed)
  s3 <- p3$run_begin(pars, save_history = TRUE)
  for (i in seq_len(nrow(data1))) {
    s3$step(i)
  }

  expect_identical(c(s1$log_likelihood, s2$log_likelihood),
                   s3$log_likelihood)
  expect_identical(s3$history$value[, , 1, ], s1$history$value)
  expect_identical(s3$history$value[, , 2, ], s2$history$value)
  expect_identical(s3$history$order[, 1, ], s1$history$order)
  expect_identical(s3$history$order[, 2, ], s2$history$order)
  expect_identical(s3$history$index, s1$history$index)
})


test_that("Can offset the initial likelihood", {
  dat <- example_sir()
  n_particles <- 42

  constant_ll <- function(pars) {
    10
  }

  set.seed(1)
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index, seed = 1L)
  ll1 <- p1$run(save_history = TRUE)

  set.seed(1)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            constant_log_likelihood = constant_ll,
                            index = dat$index, seed = 1L)
  ll2 <- p2$run(save_history = TRUE)
  expect_equal(ll2, ll1 + 10)
  expect_identical(p1$history(), p2$history())
  expect_identical(p1$state(), p2$state())
})


test_that("can save history - nested", {
  dat <- example_sir_shared()
  n_particles <- 42
  set.seed(1)

  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.2))
  constant_log_likelihood <- function(p) {
    -p$beta * 10 - p$gamma
  }

  set.seed(1)
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1)
  ll1 <- p1$run(pars)

  set.seed(1)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1,
                            constant_log_likelihood = constant_log_likelihood)
  ll2 <- p2$run(pars)
  expect_equal(ll2, ll1 - c(2.1, 3.2))
})
