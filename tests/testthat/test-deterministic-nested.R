test_that("Confirm nested particle deterministic is correct", {
  dat <- example_sir_shared()

  ## The usual compare, but add a fixed amount of noise
  compare <- function(state, observed, pars = NULL) {
    if (is.na(observed$incidence)) {
      return(NULL)
    }
    incidence_modelled <- state[1, , drop = TRUE]
    incidence_observed <- observed$incidence
    lambda <- incidence_modelled + 1e-7
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
  }

  pars <- list(list(beta = 0.2, gamma = 0.1),
               list(beta = 0.3, gamma = 0.1))

  data1 <- particle_filter_data(dat$data_raw[dat$data_raw$populations == "a", ],
                                time = "day", rate = 4)
  data2 <- particle_filter_data(dat$data_raw[dat$data_raw$populations == "b", ],
                                time = "day", rate = 4)

  p1 <- particle_deterministic$new(data1, dat$model, compare,
                                   index = dat$index)
  p2 <- particle_deterministic$new(data2, dat$model, compare,
                                   index = dat$index)
  ll1 <- p1$run(pars[[1]], save_history = TRUE, save_restart = 70)
  ll2 <- p2$run(pars[[2]], save_history = TRUE, save_restart = 70)

  p3 <- particle_deterministic$new(dat$data, dat$model, compare,
                                   index = dat$index)
  ll3 <- p3$run(pars, save_history = TRUE, save_restart = 70)

  expect_identical(ll3, c(ll1, ll2))
  expect_identical(
    array_drop(p3$history()[, , 1, , drop = FALSE], 3),
    p1$history())
  expect_identical(
    array_drop(p3$history()[, , 2, , drop = FALSE], 3),
    p2$history())
  expect_identical(
    array_drop(p3$restart_state()[, , 1, , drop = FALSE], 3),
    p1$restart_state())
  expect_identical(
    array_drop(p3$restart_state()[, , 2, , drop = FALSE], 3),
    p2$restart_state())
})


test_that("error on different population indices", {
  dat <- example_sir_shared()
  index <- function(info) {
    list(run = info$pars$beta * 10)
  }
  p <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                  index = index)
  pars <- list(
    list(beta = 0.2, gamma = 0.1),
    list(beta = 0.3, gamma = 0.1))
  expect_error(p$run(pars, save_history = TRUE),
               "index must be identical across populations")
})


test_that("Can run an mcmc on a nested model", {
  dat <- example_sir_shared()

  p <- particle_deterministic$new(dat$data, dat$model, dat$compare,
                                  dat$index)

  control <- pmcmc_control(30, save_state = TRUE, save_trajectories = TRUE,
                           save_restart = 20, rerun_every = 10)
  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = matrix(0.00026),
    proposal_varied = matrix(0.00057))

  ## Extremely basic, just test we run for now
  res <- pmcmc(pars, p, control = control)
  expect_equal(dim(res$pars), c(30, 2, 2))
  expect_equal(dim(res$probabilities), c(30, 3, 2))
  expect_equal(dim(res$state), c(5, 2, 30))
  expect_equal(dim(res$trajectories$state), c(3, 2, 30, 101))
  expect_equal(dim(res$restart$state), c(5, 2, 30, 1))
})
