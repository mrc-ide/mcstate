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
