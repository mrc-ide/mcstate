test_that("A trivial multistage filter is identical to single stage", {
  dat <- example_sir()
  ## TODO: move this into the example, will break some tests though?
  index <- function(info) {
    list(run = 5L, state = c(S = 1, I = 2, R = 3))
  }

  epochs <- list()
  pars <- dat$pars$model(dat$pars$initial())
  set.seed(1)
  filter1 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = index, seed = 1L)
  ll1 <- filter1$run(pars, save_history = TRUE)

  set.seed(1)
  filter2 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = index, seed = 1L)
  ll2 <- filter2$run_staged(pars, epochs, save_history = TRUE)

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

  epochs <- list(filter_epoch(40), filter_epoch(80))
  pars <- dat$pars$model(dat$pars$initial())
  set.seed(1)
  filter1 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = index, seed = 1L)
  ll1 <- filter1$run(pars, save_history = TRUE)

  set.seed(1)
  filter2 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = index, seed = 1L)
  ll2 <- filter2$run_staged(pars, epochs, save_history = TRUE)

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
    filter_epoch(40, state = transform_state),
    filter_epoch(80, state = transform_state))

  pars <- dat$pars$model(dat$pars$initial())
  set.seed(1)
  filter1 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = dat$index, seed = 1L)
  ll1 <- filter1$run(pars)

  set.seed(1)
  filter2 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = dat$index, seed = 1L)
  ll2 <- filter2$run_staged(pars, epochs)

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
  pars <- list(len = 10, sd = 1)
  epochs <- list(
    filter_epoch(40,
                 pars = list(len = 20, sd = 1),
                 state = transform_state),
    filter_epoch(100,
                 pars = list(len = 5, sd = 1),
                 state = transform_state))

  filter <- particle_filter$new(data, model, 42,
                                compare = compare, index = index,
                                seed = 1L)
  ## Here we just check that we can run this at all.
  expect_silent(filter$run_staged(pars, epochs, save_history = TRUE))
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
  pars <- list(len = 10)
  epochs <- list(
    filter_epoch(40,
                 pars = list(len = 20),
                 state = transform_state),
    filter_epoch(100,
                 pars = list(len = 15),
                 state = transform_state))

  filter <- new_filter()
  ll_staged <- filter$run_staged(pars, epochs, save_history = TRUE)
  h_staged <- filter$history()
  expect_identical(ll_staged, ll_20)

  ## OK, we are off here, and quite impressively so.

  ## common states are easy to check:
  expect_equal(h_staged[1:5, , ], h_20[1:5, , ])

  h_cmp <- h_20
  h_cmp[6:10, ,  1:11] <- NA
  h_cmp[9:10, , 27:51] <- NA
  expect_identical(is.na(h_staged), is.na(h_cmp))
  expect_identical(h_staged, h_cmp)
})
