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
  testthat::skip_if_not_installed("odin.dust")
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
                 transform_state = dat$transform_state),
    multistage_epoch(25,
                 pars = list(len = 15),
                 transform_state = dat$transform_state))
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
  testthat::skip_if_not_installed("odin.dust")
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
                 transform_state = dat$transform_state),
    multistage_epoch(25,
                 pars = list(len = 5, sd = 1),
                 transform_state = dat$transform_state))
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


test_that("Confirm deterministic nested multistage is correct", {
  testthat::skip_if_not_installed("odin.dust")
  dat <- example_variable()

  ## We need some multipopulation data here:
  data_raw <- data.frame(t = rep(1:50, 2),
                         observed = rnorm(100),
                         population = factor(rep(c("a", "b"), each = 50)))
  data <- particle_filter_data(data_raw, population = "population",
                               time = "t", rate = 4)
  new_filter <- function() {
    set.seed(1)
    particle_deterministic$new(data, dat$model, compare = dat$compare,
                               index = dat$index)
  }

  ## First show that the likelihood and associated history does not vary
  ## with the number of states for the simple filter; we need this
  ## property to hold to do the next tests.
  filter_10 <- new_filter()
  p_10 <- list(list(len = 10, sd = 1), list(len = 10, sd = 2))
  ll_10 <- filter_10$run(p_10, save_history = TRUE)
  h_10 <- filter_10$history()

  filter_15 <- new_filter()
  p_15 <- list(list(len = 15, sd = 1), list(len = 15, sd = 2))
  ll_15 <- filter_15$run(p_15, save_history = TRUE)
  h_15 <- filter_15$history()

  filter_20 <- new_filter()
  p_20 <- list(list(len = 20, sd = 1), list(len = 20, sd = 2))
  ll_20 <- filter_20$run(p_20, save_history = TRUE)
  h_20 <- filter_20$history()

  expect_identical(ll_10, ll_20)
  expect_identical(ll_15, ll_20)
  expect_identical(h_10, h_20[1:5, , , , drop = FALSE])
  expect_identical(h_15, h_20[1:8, , , , drop = FALSE])

  ## Now, we can set up some runs where we fiddle with the size of the
  ## model over time, and we should see that it matches the histories
  ## above:
  pars_base1 <- list(len = 10, sd = 1)
  epochs1 <- list(
    multistage_epoch(10,
                 pars = list(len = 20, sd = 1),
                 transform_state = dat$transform_state),
    multistage_epoch(25,
                 pars = list(len = 15, sd = 1),
                 transform_state = dat$transform_state))
  pars1 <- multistage_parameters(pars_base1, epochs1)

  pars_base2 <- list(len = 10, sd = 2)
  epochs2 <- list(
    multistage_epoch(10,
                 pars = list(len = 20, sd = 2),
                 transform_state = dat$transform_state),
    multistage_epoch(25,
                 pars = list(len = 15, sd = 2),
                 transform_state = dat$transform_state))
  pars2 <- multistage_parameters(pars_base2, epochs2)

  pars <- list(pars1, pars2)

  filter <- new_filter()
  ll_staged <- filter$run(pars, save_history = TRUE)
  h_staged <- filter$history()
  expect_equal(ll_staged, ll_20, tolerance = 1e-12)

  ## Common states are easy to check:
  expect_identical(h_staged[1:5, , , ], h_20[1:5, , , ])

  ## The full set requires a bit more work:
  h_cmp <- h_20
  h_cmp[6:10, , ,  1:11] <- NA
  h_cmp[9:10, , , 27:51] <- NA
  expect_identical(is.na(h_staged), is.na(h_cmp))
  expect_identical(h_staged, h_cmp)
})


test_that("Can run multistage with compiled", {
  testthat::skip_if_not_installed("odin.dust")
  dat <- example_variable()

  ## We need some multipopulation data here:
  data_raw <- data.frame(t = rep(1:50, 2),
                         observed = rnorm(100),
                         population = factor(rep(c("a", "b"), each = 50)))
  data <- particle_filter_data(data_raw, population = "population",
                               time = "t", rate = 4)

  p1 <- particle_deterministic$new(data, dat$model, compare = dat$compare,
                                   index = dat$index)
  p2 <- particle_deterministic$new(data, dat$model, compare = NULL,
                                   index = dat$index)

  ## Now, we can set up some runs where we fiddle with the size of the
  ## model over time, and we should see that it matches the histories
  ## above:
  pars_base1 <- list(len = 10, sd = 1)
  epochs1 <- list(
    multistage_epoch(10,
                 pars = list(len = 20, sd = 1),
                 transform_state = dat$transform_state),
    multistage_epoch(25,
                 pars = list(len = 15, sd = 1),
                 transform_state = dat$transform_state))
  pars1 <- multistage_parameters(pars_base1, epochs1)

  pars_base2 <- list(len = 10, sd = 2)
  epochs2 <- list(
    multistage_epoch(10,
                 pars = list(len = 20, sd = 2),
                 transform_state = dat$transform_state),
    multistage_epoch(25,
                 pars = list(len = 15, sd = 2),
                 transform_state = dat$transform_state))
  pars2 <- multistage_parameters(pars_base2, epochs2)

  pars <- list(pars1, pars2)

  ll1 <- p1$run(pars, save_history = TRUE)
  h1 <- p1$history()

  ll2 <- p2$run(pars, save_history = TRUE)
  h2 <- p2$history()

  expect_equal(ll1, ll2)
  expect_equal(h1, h2)
})


test_that("multistage model does not reuse data after initialisation", {
  testthat::skip_if_not_installed("odin.dust")
  dat <- example_variable()
  pars_base <- list(len = 10)
  epochs <- list(
    multistage_epoch(10,
                 pars = list(len = 20),
                 transform_state = dat$transform_state),
    multistage_epoch(25,
                 pars = list(len = 15),
                 transform_state = dat$transform_state))
  pars <- multistage_parameters(pars_base, epochs)

  filter <- particle_deterministic$new(dat$data, dat$model, compare = NULL,
                                       index = dat$index)
  ## This fully initialises the model, setting data:
  filter$run(pars)

  ## Set these to things that otherwise will error on use (here this
  ## will fail at model$set_data(data_split) which is what we want to
  ## avoid.  We can't achieve this via mocking because of the way that
  ## R6 organises classes and how mockery replace things
  ## https://github.com/r-lib/mockery/issues/21
  private <- r6_private(filter)
  private$data_split <- 1
  expect_silent(filter$run(pars))
})
