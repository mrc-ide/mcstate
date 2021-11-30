test_that("A trivial multistage filter is identical to single stage", {
  dat <- example_sir()

  epochs <- list()
  pars <- dat$pars$model(dat$pars$initial())
  set.seed(1)
  filter1 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = dat$index, seed = 1L)
  ll1 <- filter1$run(pars)

  set.seed(1)
  filter2 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = dat$index, seed = 1L)
  ll2 <- filter2$run_staged(pars, epochs)

  expect_identical(ll1, ll2)
  expect_identical(
    r6_private(filter1)$last_model$rng_state(),
    r6_private(filter2)$last_model$rng_state())
})


test_that("An effectless multistage filter is identical to single stage", {
  dat <- example_sir()

  epochs <- list(filter_epoch(40), filter_epoch(80))
  pars <- dat$pars$model(dat$pars$initial())
  set.seed(1)
  filter1 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = dat$index, seed = 1L)
  ll1 <- filter1$run(pars)

  set.seed(1)
  filter2 <- particle_filter$new(dat$data, dat$model, 42, dat$compare,
                                 index = dat$index, seed = 1L)
  ll2 <- filter2$run_staged(pars, epochs)

  expect_identical(ll1, ll2)
  expect_identical(
    r6_private(filter1)$last_model$rng_state(),
    r6_private(filter2)$last_model$rng_state())
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

  transform_trajectories <- function(y, model) {
    y[1:5, ]
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
                 state = transform_state,
                 trajectories = transform_trajectories)
    filter_epoch(100,
                 pars = list(len = 5, sd = 1),
                 state = transform_state,
                 trajectories = transform_trajectories))

  filter <- particle_filter$new(data, model, 42, compare = compare, seed = 1L)
  filter$run_staged(pars, epochs)

  ## What can I usefully test for here? nothing yet.
})
