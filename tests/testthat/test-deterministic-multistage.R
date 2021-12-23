test_that("Trivial multistage deterministic same as simple", {
  dat <- example_sir()
  index <- function(info) {
    list(run = 5L, state = c(S = 1, I = 2, R = 3))
  }

  pars_base <- list(beta = 0.2, gamma = 0.1, compare = list(exp_noise = Inf))
  pars <- multistage_parameters(pars_base, list())

  set.seed(1)
  filter1 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                        index = index)
  ll1 <- filter1$run(pars_base, save_history = TRUE)

  set.seed(2)
  filter2 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                        index = index)
  ll2 <- filter2$run(pars, save_history = TRUE)

  expect_identical(ll1, ll2)

  ## expect_identical(
  ##   r6_private(filter1)$last_model$rng_state(),
  ##   r6_private(filter2)$last_model$rng_state())

  ## Not totally sure how we don't have rownames here.
  expect_identical(filter1$history(), filter2$history())
})


test_that("An effectless multistage filter is identical to single stage", {
  dat <- example_sir()
  index <- function(info) {
    list(run = 5L, state = c(S = 1, I = 2, R = 3))
  }

  epochs <- list(multistage_epoch(10), multistage_epoch(20))
  pars_base <- list(beta = 0.2, gamma = 0.1, compare = list(exp_noise = Inf))
  pars <- multistage_parameters(pars_base, epochs)
  restart <- c(5, 15, 25)

  set.seed(1)
  filter1 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                        index = index)
  ll1 <- filter1$run(pars_base, save_history = TRUE, save_restart = restart)

  set.seed(1)
  filter2 <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                        index = index)
  ll2 <- filter2$run(pars, save_history = TRUE, save_restart = restart)

  ## To get these to be identical we must loop over the data in the
  ## other order (by pars-within-time, then time); see
  ## test-deterministic.R for a similar problem (but note that the
  ## actual states here are spot on)
  expect_equal(ll1, ll2, tolerance = 1e-11)
  ## expect_identical(
  ##   r6_private(filter1)$last_model$rng_state(),
  ##   r6_private(filter2)$last_model$rng_state())

  expect_identical(filter1$history(), filter2$history())
  expect_identical(filter1$restart_state(), filter2$restart_state())
})


test_that("multistage, dimension changing, model agrees with single stage", {
  dat <- example_variable()
  new_filter <- function() {
    set.seed(1)
    particle_deterministic$new(dat$data, dat$model, compare = dat$compare,
                               index = dat$index)
  }

  ## First show that the likelihood and associated history does not vary
  ## with the number of states for the simple filter; we need this
  ## property to hold to do the next tests.
  filter_10 <- new_filter()
  ll_10 <- filter_10$run(list(len = 10), save_history = TRUE)
  h_10 <- filter_10$history()

  filter_15 <- new_filter()
  ll_15 <- filter_15$run(list(len = 15), save_history = TRUE)
  h_15 <- filter_15$history()

  filter_20 <- new_filter()
  ll_20 <- filter_20$run(list(len = 20), save_history = TRUE)
  h_20 <- filter_20$history()

  expect_identical(ll_10, ll_20)
  expect_identical(ll_15, ll_20)
  expect_identical(h_10, h_20[1:5, , , drop = FALSE])
  expect_identical(h_15, h_20[1:8, , , drop = FALSE])

  ## Now, we can set up some runs where we fiddle with the size of the
  ## model over time, and we should see that it matches the histories
  ## above:
  pars_base <- list(len = 10)
  epochs <- list(
    multistage_epoch(10,
                 pars = list(len = 20),
                 transform_state = dat$transform_state_deterministic),
    multistage_epoch(25,
                 pars = list(len = 15),
                 transform_state = dat$transform_state_deterministic))
  pars <- multistage_parameters(pars_base, epochs)

  ## Looks like I have a slight issue with creating the new model
  ## state, something to fix up next, but probably not insurmountable.

  ## Most of the issues here are that we set this up to work for any
  ## number of parameters at once, and that is just not working well
  ## here.
  filter <- new_filter()
  ll_staged <- filter$run(pars, save_history = TRUE)
  h_staged <- filter$history()
  expect_equal(ll_staged, ll_20, tolerance = 1e-10)

  ## Common states are easy to check:
  expect_equal(h_staged[1:5, , ], h_20[1:5, , ])

  ## The full set requires a bit more work:
  h_cmp <- h_20
  h_cmp[6:10, ,  1:11] <- NA
  h_cmp[9:10, , 27:51] <- NA
  expect_identical(is.na(h_staged), is.na(h_cmp))
  expect_identical(h_staged, h_cmp)
})


test_that("Can't save restart when size changes", {
  dat <- example_variable()

  ## There's going to be some work here to update this so that it's
  ## easy to work with.  Possibly we just set a blank epoch at the
  ## beginning - that's not terrible and is at least symmetrical.
  ##
  ## There's some pretty major work in pmcmc to get this sorted though
  ## as we need to generate all of this out of the mcmc parameters,
  ## and that's its own challenge.
  pars_base <- list(len = 10, sd = 1)
  epochs <- list(
    multistage_epoch(10,
                 pars = list(len = 20, sd = 1),
                 transform_state = dat$transform_state_deterministic),
    multistage_epoch(25,
                 pars = list(len = 5, sd = 1),
                 transform_state = dat$transform_state_deterministic))
  pars <- multistage_parameters(pars_base, epochs)

  filter <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                       index = dat$index)

  ## Fine getting restart from the last stage
  filter$run(pars, save_restart = c(30, 40, 50))
  expect_equal(dim(filter$restart_state()), c(5, 1, 3))

  ## Or some other stage
  filter$run(pars, save_restart = c(11, 15, 20))
  expect_equal(dim(filter$restart_state()), c(20, 1, 3))

  ## But can't do anything here:
  err <- expect_error(
    filter$run(pars, save_restart = c(5, 15, 30)),
    "Restart state varies in size over the simulation")
  expect_match(err$message, "2: 20 rows", all = FALSE)
})


test_that("Parameters handling works", {
  expect_error(
    particle_deterministic_pars(list()),
    "At least one parameter set required")
  expect_equal(
    particle_deterministic_pars(list(list(a = 1))),
    list(list(a = 1)))
  expect_equal(
    particle_deterministic_pars(list(list(a = 1), list(a = 2))),
    list(list(a = 1), list(a = 2)))

  ## To test how things combine we need two of each types

  f <- function(...) identity(..1)

  ## Basic parameters
  p1a <- list(a = 1)
  p1b <- list(a = 10)

  ## Trivial multistage parameters
  p2a <- multistage_parameters(list(a = 2), list())
  p2b <- multistage_parameters(list(a = 20), list())

  ## With a single epoch changing parameters, but not transform
  p3a <- multistage_parameters(list(a = 3),
                               list(multistage_epoch(10, p1a)))
  p3b <- multistage_parameters(list(a = 30),
                               list(multistage_epoch(10, p1b)))

  ## With a single epoch chnging both parameters and transform
  p4a <- multistage_parameters(list(a = 4),
                               list(multistage_epoch(10, p1a, f)))
  p4b <- multistage_parameters(list(a = 40),
                               list(multistage_epoch(10, p1b, f)))

  ## With a single epoch chnging transform, but not parameters
  p5a <- multistage_parameters(list(a = 5),
                               list(multistage_epoch(10, NULL, f)))
  p5b <- multistage_parameters(list(a = 60),
                               list(multistage_epoch(10, NULL, f)))

  others <- c("start", "transform_state")

  res <- particle_deterministic_pars(list(p2a))
  expect_s3_class(res, "multistage_parameters")
  expect_length(res, 1)
  expect_equal(res[[1]]$pars[[1]], p2a[[1]]$pars)
  expect_equal(res[[1]][others], p2a[[1]][others])

  res <- particle_deterministic_pars(list(p2a, p2b))
  expect_s3_class(res, "multistage_parameters")
  expect_length(res, 1) # still one stage
  expect_equal(res[[1]]$pars[[1]], p2a[[1]]$pars)
  expect_equal(res[[1]]$pars[[2]], p2b[[1]]$pars)
  expect_equal(res[[1]][others], p2a[[1]][others])

  res <- particle_deterministic_pars(list(p3a, p3b))
  expect_s3_class(res, "multistage_parameters")
  expect_length(res, 2)
  expect_equal(lapply(res, function(x) x$pars),
               Map(function(a, b) list(a$pars, b$pars), p3a, p3b))
  expect_equal(lapply(res, function(x) x[others]),
               lapply(p3a, function(x) x[others]))

  res <- particle_deterministic_pars(list(p4a, p4b))
  expect_s3_class(res, "multistage_parameters")
  expect_length(res, 2)
  expect_equal(lapply(res, function(x) x$pars),
               Map(function(a, b) list(a$pars, b$pars), p4a, p4b))
  expect_equal(lapply(res, function(x) x[others]),
               lapply(p4a, function(x) x[others]))

  res <- particle_deterministic_pars(list(p5a, p5b))
  expect_s3_class(res, "multistage_parameters")
  expect_length(res, 2)
  expect_equal(lapply(res, function(x) x$pars),
               list(list(p5a[[1]]$pars, p5b[[1]]$pars), NULL))
  expect_equal(lapply(res, function(x) x[others]),
               lapply(p5a, function(x) x[others]))

  expect_error(
    particle_deterministic_pars(list(p1a, p2a)),
    "'pars' must be either all multistage or all non-multistage")
  expect_error(
    particle_deterministic_pars(list(p2a, p3a)),
    "Incompatible numbers of stages in pars: found 1, 2 stages")

  p_start <- multistage_parameters(list(a = 5), list(multistage_epoch(11)))
  expect_error(
    particle_deterministic_pars(list(p3a, p_start)),
    "Incompatible 'start' time at phase 2")

  expect_error(
    particle_deterministic_pars(list(p3a, p4a)),
    "Incompatible 'transform_state' at phase 2")

  p_no_pars <- multistage_parameters(list(a = 5), list(multistage_epoch(10)))
  expect_error(
    particle_deterministic_pars(list(p3a, p_no_pars)),
    "Incompatible 'pars' at phase 2")
})
