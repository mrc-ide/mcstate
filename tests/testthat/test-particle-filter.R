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
    list(state = c(1000, pars$I0, 0, 0, 0),
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
    list(state = rbind(1000, i0, 0, 0, 0, deparse.level = 0),
         step = step)
  }
  compare <- function(...) {
    NULL
  }

  p <- particle_filter$new(data, dat$model, n_particles, compare,
                           index = dat$index, initial = initial, seed = 1L)
  pars <- list(I0 = 10, step_mean = 20, step_offset = 400)
  set.seed(1)
  ll <- p$run(pars, save_history = TRUE)
  expect_equal(ll, 0)

  mod <- p$model$new(list(), step = 0, n_particles = n_particles, seed = 1L)
  tmp <- initial(NULL, n_particles, pars)
  mod$set_state(tmp$state, tmp$step)
  expect_equal(p$history()[, , 1], mod$state()[1:3, ])
  expect_equal(mod$step(), max(tmp$step))
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


test_that("prevent using compiled compare where model does not support it", {
  dat <- example_sir()
  n_particles <- 100
  model <- dust::dust_example("walk")
  expect_error(
    particle_filter$new(dat$data, model, n_particles, NULL),
    "Your model does not have a built-in 'compare' function")
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
  res <- obj$fork(list(beta = 0.15))
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

test_that("particle filter state nested - errors", {
  dat <- example_sir_shared()
  n_particles <- 42
  set.seed(1)

  pars <- list(list(beta = 0.2, gamma = 0.1),
                               list(beta = 0.3, gamma = 0.1))

  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)

  expect_error(p$run(pars[1]), "the length")
})

test_that("nested pf with initial", {
  initial <- function(info, n_particles, pars) {
    list(step = pars[[2]]$initial)
  }

  dat <- example_sir_shared()
  n_particles <- 42
  data <- dat$data
  offset <- 400
  data[c("step_start", "step_end")] <-
    data[c("step_start", "step_end")] + offset
  data$step_start[c(1, 101)] <- 0

  pars <- list(list(beta = 0.2, gamma = 0.1),
                               list(beta = 0.3, gamma = 0.1))

  ## The usual version:
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index)
  set.seed(1)
  ll1 <- p1$run(pars)

  ## Tuning the start date
  p2 <- particle_filter$new(data, dat$model, n_particles, dat$compare,
                            index = dat$index, initial = initial)
  set.seed(1)
  pars <- list(list(beta = 0.2, gamma = 0.1, initial = as.integer(offset)),
               list(beta = 0.3, gamma = 0.1, initial = as.integer(offset)))
  ll2 <- p2$run(pars)
  expect_identical(ll1, ll2)

  ## Running from the beginning is much worse:
  set.seed(1)
  ll3 <- p2$run(list(pars, list(initial = 0L)))
  expect_true(all(ll3 < ll1))
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
  res <- obj$fork(pars)
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