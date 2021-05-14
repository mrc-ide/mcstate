context("if2")

test_that("Can run IF2", {
  pars <- if2_parameters$new(
            list(if2_parameter("beta", 0.15, min = 0, max = 1),
                 if2_parameter("gamma", 0.05, min = 0, max = 1)))
  pars_sd <- c(beta = 0.02, gamma = 0.02) # not a list!

  dat <- example_sir()

  filter <- if2$new(pars, dat$data, dat$model, dat$compare, NULL, dat$index)

  iterations <- 50
  cooling_target <- 0.5
  n_par_sets <- 100
  set.seed(1)

  filter$run(pars_sd, iterations, n_par_sets, cooling_target)

  expect_equal(length(filter$log_likelihood()), iterations)
  expect_equal(dim(filter$pars_series()),
               c(length(pars$names()), n_par_sets, iterations))

  filter$plot()
  filter$plot("beta")
  filter$plot("gamma")

  n_particles <- 100
  filter$sample(n_particles)
})

