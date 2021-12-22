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
