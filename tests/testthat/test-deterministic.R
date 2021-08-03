context("deterministic")

test_that("Can run the deterministic filter", {
  dat <- example_sir()
  p <- particle_nofilter$new(dat$data, dat$model, dat$compare, dat$index)
  ll <- p$run(dat$pars$model(dat$pars$initial()))
})


test_that("Can control starting point of simulation", {
  initial <- function(info, n_particles, pars) {
    list(step = pars$initial)
  }

  dat <- example_sir()

  data <- dat$data
  offset <- 400
  data[c("step_start", "step_end")] <-
    data[c("step_start", "step_end")] + offset
  data$step_start[[1]] <- 0

  pars <- dat$pars$model(dat$pars$initial())

  ## The usual version:
  p1 <- particle_nofilter$new(dat$data, dat$model, dat$compare,
                              index = dat$index)
  set.seed(1)
  ll1 <- p1$run(pars)

  ## Tuning the start date
  p2 <- particle_nofilter$new(data, dat$model, dat$compare,
                              index = dat$index, initial = initial)
  set.seed(1)
  ll2 <- p2$run(c(pars, list(initial = as.integer(offset))))
  expect_identical(ll1, ll2)
})
