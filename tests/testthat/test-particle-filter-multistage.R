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
  ## TODO: move this into the example, will break some tests though?
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

  ## TODO: see above
  index <- function(info) {
    list(run = 5L, state = c(S = 1, I = 2, R = 3))
  }

  epochs <- list(multistage_epoch(10), multistage_epoch(20))
  pars_base <- dat$pars$model(dat$pars$initial())
  pars <- multistage_parameters(pars_base, epochs)
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


test_that("Can transform state in the model", {
  dat <- example_sir()

  ## Something that will leave an effect; if we add new individuals at
  ## these time points then we'll see the population size increase.
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

  set.seed(1)
  filter2 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = dat$index, seed = 1L)
  ll2 <- filter2$run(pars)

  expect_false(ll1 == ll2)
  expect_equal(colSums(filter1$state(1:3)), rep(1010, 42))
  expect_equal(colSums(filter2$state(1:3)), rep(1050, 42))
})


test_that("Can transform state size", {
  model <- dust::dust_example("variable")
  data <- particle_filter_data(data.frame(time = 1:50, observed = rnorm(50)),
                               "time", 4)
  ## Nonsense model
  compare <- function(state, observed, pars) {
    apply(dnorm(state, log = TRUE), 2, max)
  }

  index <- function(info) {
    i <- seq(1, info$len, by = 2L)
    names(i) <- letters[i]
    list(run = i, state = i)
  }

  transform_state <- function(y, model_old, model_new) {
    n_old <- model_old$pars()$len
    n_new <- model_new$pars()$len
    if (n_new > n_old) {
      y <- rbind(y, matrix(0, n_new - n_old, ncol(y)))
    } else {
      y <- y[seq_len(n_new), ]
    }
    y
  }

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
                 transform_state = transform_state),
    multistage_epoch(25,
                 pars = list(len = 5, sd = 1),
                 transform_state = transform_state))
  pars <- multistage_parameters(pars_base, epochs)

  filter <- particle_filter$new(data, model, 42,
                                compare = compare, index = index,
                                seed = 1L)
  ## Here we just check that we can run this at all.
  expect_silent(filter$run(pars, save_history = TRUE))
})


test_that("multistage, dimension changing, model agrees with single stage", {
  ## A small, very silly, model designed to help work with the
  ## multistage filter.  We have a model we can change the dimensions of
  ## without changing the way that the random number draws will work
  ## because only the first entry will be stochastic.
  model <- odin.dust::odin_dust({
    len <- user(integer = TRUE)
    update(x[1]) <- x[1] + rnorm(0, 0.1)
    update(x[2:len]) <- i + step / 10
    initial(x[]) <- 0
    dim(x) <- len
  }, verbose = FALSE)

  data <- particle_filter_data(data.frame(time = 1:50, observed = rnorm(50)),
                               "time", 4)
  ## Nonsense model
  compare <- function(state, observed, pars) {
    dnorm(state - observed$observed, log = TRUE)
  }

  index <- function(info) {
    i <- seq(1, info$len, by = 2L)
    names(i) <- letters[i]
    list(run = 1L, state = i)
  }

  transform_state <- function(y, model_old, model_new) {
    n_old <- model_old$pars()$len
    n_new <- model_new$pars()$len
    if (n_new > n_old) {
      y <- rbind(y, matrix(0, n_new - n_old, ncol(y)))
    } else {
      y <- y[seq_len(n_new), ]
    }
    y
  }

  new_filter <- function() {
    set.seed(1)
    particle_filter$new(data, model, 42,
                        compare = compare, index = index,
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
                 transform_state = transform_state),
    multistage_epoch(25,
                 pars = list(len = 15),
                 transform_state = transform_state))
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
  dat <- example_sir()

  ## TODO: see above
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
