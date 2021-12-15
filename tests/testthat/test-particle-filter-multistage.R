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
    r6_private(filter1)$last_model$rng_state(),
    r6_private(filter2)$last_model$rng_state())

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
    r6_private(filter1)$last_model$rng_state(),
    r6_private(filter2)$last_model$rng_state())

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
  ## Here we just check that we can run this at all.
  expect_silent(filter$run(pars, save_history = TRUE))
})


test_that("multistage, dimension changing, model agrees with single stage", {
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


test_that("All times must be found in the data", {
  skip("now obsolete")
  dat <- example_sir()
  index <- function(info) {
    list(run = 5L, state = c(S = 1, I = 2, R = 3))
  }

  pars_base <- dat$pars$model(dat$pars$initial())
  epochs <- list(multistage_epoch(10.5),
                 multistage_epoch(20),
                 multistage_epoch(200))
  filter <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                index = index, seed = 1L)
  expect_error(
    filter$run(multistage_parameters(pars_base, epochs[1:2])),
    "Could not map epoch to filter time: error for 1")
  expect_error(
    filter$run(multistage_parameters(pars_base, epochs)),
    "Could not map epoch to filter time: error for 1, 3")
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


test_that("Can run a multistage filter from part way through", {
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

  t_min <- 35
  step_min <- t_min * attr(dat$data, "rate")
  data <- dat$data
  data <- dat$data[dat$data$time_end > t_min, ]

  initial <- function(info, n_particles, pars) {
    list(state = rep(0, info$len), step = step_min)
  }

  filter <- particle_filter$new(data, dat$model, 42,
                                compare = dat$compare, index = dat$index,
                                initial = initial, seed = 1L)
  filter$run(pars, save_history = TRUE)

  ## Here we just check that we can run this at all.
  expect_silent(filter$run(pars, save_history = TRUE))
})
