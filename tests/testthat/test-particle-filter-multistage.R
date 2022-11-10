test_that("can construct multistage_epoch objects", {
  expect_equal(
    multistage_epoch(10),
    structure(list(start = 10, pars = NULL,
                   transform_state = transform_state_identity),
              class = "multistage_epoch"))
  transform <- function(x, ...) {
    x + 1
  }
  expect_equal(
    multistage_epoch(10, transform_state = transform),
    structure(list(start = 10, pars = NULL, transform_state = transform),
              class = "multistage_epoch"))
  expect_equal(
    multistage_epoch(10, list(a = 1), transform),
    structure(list(start = 10, pars = list(a = 1), transform_state = transform),
              class = "multistage_epoch"))

  expect_error(
    multistage_epoch(10, list(a = 1), TRUE),
    "'transform_state' must be a function")
})


test_that("Can construct trivial multistage_parameters object", {
  base <- list(a = 1)
  res <- multistage_parameters(base, list())
  expect_equal(
    res,
    structure(
      list(list(pars = base, start = NULL, transform_state = NULL)),
      class = "multistage_parameters"))
})


test_that("Can construct multistage_parameters object with 2 stages", {
  base <- list(a = 1)
  epochs <- list(
    multistage_epoch(20, list(a = 2)),
    multistage_epoch(30, list(a = 3)))

  res <- multistage_parameters(base, epochs)
  expect_length(res, 3)
  expect_s3_class(res, "multistage_parameters")

  expect_equal(res[[1]],
               list(pars = base, start = NULL, transform_state = NULL))
  expect_equal(res[[2]],
               list(pars = list(a = 2),
                    start = 20,
                    transform_state = transform_state_identity))
  expect_equal(res[[3]],
               list(pars = list(a = 3),
                    start = 30,
                    transform_state = transform_state_identity))
})


test_that("parameter changes must be consistent", {
  base <- list(a = 1)
  epochs <- list(
    multistage_epoch(20),
    multistage_epoch(30, list(a = 3)))
  expect_error(
    multistage_parameters(base, epochs),
    "If changing parameters, all epochs must contain 'pars'")
})


test_that("all epoch entries must be multistage_epoch objects", {
  base <- list(a = 1)
  epochs <- list(
    multistage_epoch(20),
    NULL)
  expect_error(
    multistage_parameters(base, epochs),
    "Expected all elements of 'epochs' to be 'multistage_epoch' objects")
})


test_that("A trivial multistage filter is identical to single stage", {
  dat <- example_sir()

  ## we need a named index throughout here, and can't use the one in
  ## the example.
  index <- function(info) {
    list(run = 5L, state = c(S = 1, I = 2, R = 3))
  }
  pars_base <- dat$pars$model(dat$pars$initial())
  pars <- multistage_parameters(pars_base, list())

  set.seed(1)
  filter1 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = index, seed = 1L)
  ll1 <- filter1$run(pars_base, save_history = TRUE)

  set.seed(1)
  filter2 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = index, seed = 1L)
  ll2 <- filter2$run(pars, save_history = TRUE)

  expect_identical(ll1, ll2)
  expect_identical(
    r6_private(filter1)$last_model[[1]]$rng_state(),
    r6_private(filter2)$last_model[[1]]$rng_state())

  expect_identical(filter1$history(), filter2$history())
})


test_that("An effectless multistage filter is identical to single stage", {
  dat <- example_sir()
  index <- function(info) {
    list(run = 5L, state = c(S = 1, I = 2, R = 3))
  }

  epochs <- list(multistage_epoch(10), multistage_epoch(20))
  pars_base <- dat$pars$model(dat$pars$initial())
  pars <- multistage_parameters(pars_base, epochs)
  restart <- c(5, 15, 25)

  set.seed(1)
  filter1 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = index, seed = 1L)
  ll1 <- filter1$run(pars_base, save_history = TRUE, save_restart = restart)

  set.seed(1)
  filter2 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = index, seed = 1L)
  ll2 <- filter2$run(pars, save_history = TRUE, save_restart = restart)

  expect_identical(ll1, ll2)
  expect_identical(
    r6_private(filter1)$last_model[[1]]$rng_state(),
    r6_private(filter2)$last_model[[3]]$rng_state())

  expect_identical(filter1$history(), filter2$history())
  expect_identical(filter1$restart_state(), filter2$restart_state())
})


test_that("Can transform state in the model", {
  dat <- example_sir()

  ## Something that will leave an effect; if we add new individuals at
  ## these time points then we'll see the population size increase,
  ## both at the end of the simulation and in the restart data.
  transform_state <- function(y, model_old, model_new) {
    y[2, ] <- y[2, ] + 20
    y
  }

  epochs <- list(
    multistage_epoch(10, transform_state = transform_state),
    multistage_epoch(20, transform_state = transform_state))

  pars_base <- dat$pars$model(dat$pars$initial())
  pars <- multistage_parameters(pars_base, epochs)

  set.seed(1)
  filter1 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = dat$index, seed = 1L)
  ll1 <- filter1$run(pars_base)

  restart <- c(5, 15, 25)
  set.seed(1)
  filter2 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = dat$index, seed = 1L)
  ll2 <- filter2$run(pars, save_restart = restart)

  expect_false(ll1 == ll2)
  expect_equal(colSums(filter1$state(1:3)), rep(1010, 42))
  expect_equal(colSums(filter2$state(1:3)), rep(1050, 42))

  r <- filter2$restart_state()
  expect_equal(dim(r), c(5, 42, 3))
  expect_equal(
    apply(r[1:3, , ], 2:3, sum),
    matrix(c(1010, 1030, 1050), 42, 3, byrow = TRUE))
})


test_that("Can transform state size", {
  testthat::skip_if_not_installed("odin.dust")
  dat <- example_variable()

  pars_base <- list(len = 10, sd = 1)
  epochs <- list(
    multistage_epoch(10,
                 pars = list(len = 20, sd = 1),
                 transform_state = dat$transform_state),
    multistage_epoch(25,
                 pars = list(len = 5, sd = 1),
                 transform_state = dat$transform_state))
  pars <- multistage_parameters(pars_base, epochs)

  filter <- particle_filter$new(dat$data, dat$model, 42,
                                compare = dat$compare, index = dat$index,
                                seed = 1L)
  ## Here we just check that we can run this at all.
  expect_silent(filter$run(pars, save_history = TRUE))
})


test_that("multistage, dimension changing, model agrees with single stage", {
  testthat::skip_if_not_installed("odin.dust")
  dat <- example_variable()
  new_filter <- function() {
    set.seed(1)
    particle_filter$new(dat$data, dat$model, 42,
                        compare = dat$compare, index = dat$index,
                        seed = 1L)
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
  expect_identical(h_10, h_20[1:5, , ])
  expect_identical(h_15, h_20[1:8, , ])

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

  filter <- new_filter()
  ll_staged <- filter$run(pars, save_history = TRUE)
  h_staged <- filter$history()
  expect_identical(ll_staged, ll_20)

  ## Common states are easy to check:
  expect_equal(h_staged[1:5, , ], h_20[1:5, , ])

  ## The full set requires a bit more work:
  h_cmp <- h_20
  h_cmp[6:10, ,  1:11] <- NA
  h_cmp[9:10, , 27:51] <- NA
  expect_identical(is.na(h_staged), is.na(h_cmp))
  expect_identical(h_staged, h_cmp)
})


test_that("Require named index for history-saving multistage filter", {
  dat <- example_sir()
  index <- function(info) {
    list(run = 5L, state = 1:3)
  }

  pars_base <- dat$pars$model(dat$pars$initial())
  epochs <- list(multistage_epoch(10),
                 multistage_epoch(20))
  filter <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                index = index, seed = 1L)
  expect_error(
    filter$run(multistage_parameters(pars_base, epochs), save_history = TRUE),
    "Named index required")
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

  filter <- particle_filter$new(dat$data, dat$model, 42,
                                compare = dat$compare, index = dat$index,
                                seed = 1L)

  ## Fine getting restart from the last stage
  filter$run(pars, save_restart = c(30, 40, 50))
  expect_equal(dim(filter$restart_state()), c(5, 42, 3))

  ## Or some other stage
  filter$run(pars, save_restart = c(11, 15, 20))
  expect_equal(dim(filter$restart_state()), c(20, 42, 3))

  ## But can't do anything here:
  err <- expect_error(
    filter$run(pars, save_restart = c(5, 15, 30)),
    "Restart state varies in size over the simulation")
  expect_match(err$message, "2: 20 rows", all = FALSE)
})


test_that("Prevent restart on epoch changes", {
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

  filter <- particle_filter$new(dat$data, dat$model, 42,
                                compare = dat$compare, index = dat$index,
                                seed = 1L)
  expect_error(
    filter$run(pars, save_restart = 10),
    "save_restart cannot include epoch change: error for 10")
  expect_error(
    filter$run(pars, save_restart = c(5, 10, 15)),
    "save_restart cannot include epoch change: error for 10")
  expect_error(
    filter$run(pars, save_restart = c(5, 10, 15, 20, 25, 30)),
    "save_restart cannot include epoch change: error for 10, 25")
})


test_that("Gracefully cope with early exit", {
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

  filter <- particle_filter$new(dat$data, dat$model, 42,
                                compare = dat$compare, index = dat$index,
                                seed = 1L)

  expect_equal(
    filter$run(pars, save_restart = 40, min_log_likelihood = -20),
    -Inf)

  ## Our restart state is consistent, even if it is junk:
  expect_equal(filter$restart_state(), array(NA_real_, c(5, 42, 1)))
})


## epochs:         t2        t3
## data:
##   |----------|                    1 only
##        |-----------|              1-2, t2 must exist
##              |--------------|     1-3, t2 and t3 must exist
##                    |--|           2 only
##                        |-------|  2-3, t3 must exist
##                              |--| 3 only
test_that("Can filter multistage parameters based on data", {
  ## We need a little helper here to remove the tedium
  f <- function(t, d) {
    base <- list(i = 1)
    epochs <- lapply(seq_along(t), function(i)
      multistage_epoch(t[[i]], pars = list(i = i + 1)))
    p <- multistage_parameters(base, epochs)
    res <- filter_check_times(p, d, NULL)
    vapply(res, function(x) c(x$pars$i, x$time_index), numeric(2))
  }

  d <- particle_filter_data(data.frame(t = 11:30, value = runif(20)),
                            "t", 4)

  expect_equal(f(integer(0), d), cbind(c(1, 20)))
  ## Changes all before any data; use last
  expect_equal(f(1:3, d), cbind(c(4, 20)))
  ## Changes all after any data; use first
  expect_equal(f(41:43, d), cbind(c(1, 20)))

  ## An epoch that starts after the data ends is ignored:
  expect_equal(f(c(20, 60), d), cbind(c(1, 10), c(2, 20)))
  ## An epoch that finishes before the data starts is ignored
  expect_equal(f(c(5, 20), d), cbind(c(2, 10), c(3, 20)))

  ## The common case:
  expect_equal(f(c(15, 25), d), cbind(c(1, 5), c(2, 15), c(3, 20)))

  ## Error case; this can't really be triggered until we get non-unit
  ## time changes supported:
  expect_error(f(c(5, 20.5), d),
               "Could not map epoch to filter time: error for stage 2")
  expect_error(f(c(5, 6, 20.5, 25), d),
               "Could not map epoch to filter time: error for stage 3")

  ## Special case where things are against the first boundary:
  expect_equal(f(10, d), cbind(c(2, 20)))
  expect_equal(f(9, d), cbind(c(2, 20)))
})


## This is not a great test, but covers the key bits; that we filter
## down the parameters to the right set, in particular.
test_that("Can run a multistage filter from part way through", {
  testthat::skip_if_not_installed("odin.dust")
  dat <- example_variable()

  pars_base <- list(len = 10, sd = 1)
  epochs <- list(
    multistage_epoch(10,
                 pars = list(len = 20, sd = 1),
                 transform_state = dat$transform_state),
    multistage_epoch(25,
                 pars = list(len = 5, sd = 1),
                 transform_state = dat$transform_state))
  pars <- multistage_parameters(pars_base, epochs)

  t_min <- 35
  step_min <- t_min * attr(dat$data, "rate")
  data <- dat$data
  data <- dat$data[dat$data$t_end > t_min, ]

  initial <- function(info, n_particles, pars) {
    rep(0, info$len)
  }

  filter <- particle_filter$new(data, dat$model, 42,
                                compare = dat$compare, index = dat$index,
                                initial = initial, seed = 1L)
  filter$run(pars, save_history = TRUE)

  ## No problems stitching together the history, now all in a single
  ## case.
  h <- filter$history()
  expect_equal(dim(h), c(3, 42, nrow(data) + 1))
  expect_false(any(is.na(h)))

  p <- filter_check_times(pars, data, NULL)
  expect_length(p, 1)
  expect_equal(p[[1]]$pars$len, 5)
})



test_that("Confirm nested filter is correct", {
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
    particle_filter$new(data, dat$model, 42,
                        compare = dat$compare, index = dat$index,
                        seed = 1L)
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
  expect_identical(h_10, h_20[1:5, , , ])
  expect_identical(h_15, h_20[1:8, , , ])

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
  expect_identical(ll_staged, ll_20)

  ## Common states are easy to check:
  expect_equal(h_staged[1:5, , , ], h_20[1:5, , , ])

  ## The full set requires a bit more work:
  h_cmp <- h_20
  h_cmp[6:10, , ,  1:11] <- NA
  h_cmp[9:10, , , 27:51] <- NA
  expect_identical(is.na(h_staged), is.na(h_cmp))
  expect_identical(h_staged, h_cmp)
})


test_that("Can't change numbers of stages after creation", {
  dat <- example_sir()

  index <- function(info) {
    list(run = 5L, state = c(S = 1, I = 2, R = 3))
  }
  pars_base <- dat$pars$model(dat$pars$initial())
  pars0 <- pars_base
  pars1 <- multistage_parameters(pars_base, list())
  pars2 <- multistage_parameters(pars_base, list(multistage_epoch(10)))

  filter1 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = index, seed = 1L)
  filter1$run(pars0)
  expect_silent(filter1$run(pars1))
  expect_error(
    filter1$run(pars2),
    "Expected single-stage parameters (but given one with 2 stages)",
    fixed = TRUE)

  filter2 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = index, seed = 1L)
  filter2$run(pars2)
  expect_error(
    filter2$run(pars0),
    "Expected multistage_pars with 2 stages (but given one with 1)",
    fixed = TRUE)
  expect_error(
    filter2$run(pars1),
    "Expected multistage_pars with 2 stages (but given one with 1)",
    fixed = TRUE)
})


test_that("can run a particle filter over a subset of data, twice", {
  testthat::skip_if_not_installed("odin.dust")
  dat <- example_variable()

  index <- function(info) {
    list(run = 5L, state = c(S = 1, I = 2, R = 3))
  }

  pars_base <- list(len = 10, sd = 1)
  epochs <- list(
    multistage_epoch(10,
                 pars = list(len = 20, sd = 1),
                 transform_state = dat$transform_state),
    multistage_epoch(25,
                 pars = list(len = 5, sd = 1),
                 transform_state = dat$transform_state))
  pars <- multistage_parameters(pars_base, epochs)

  data <- subset(dat$data, t_start >= 30)
  filter <- particle_filter$new(data, dat$model, 42, dat$compare,
                                index = index, seed = 1L)
  filter$run(pars)
  ## Test of some internal bookkeeping
  expect_length(r6_private(filter)$last_model, 1)
  expect_equal(r6_private(filter)$last_stages, 3)
  ## Previously this failed; confirm we can run it a second time!
  expect_silent(filter$run(pars))
})
