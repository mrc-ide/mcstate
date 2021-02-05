context("smc2")

test_that("Can run smc2", {
  dat <- example_sir()
  n_particles <- 42
  n_parameter_sets <- 20
  f <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  p <- smc2_parameters$new(
    list(
      smc2_parameter("beta",
                     function(n) runif(n, 0, 1),
                     function(x) dunif(x, 0, 1, log = TRUE),
                     min = 0, max = 1),
      smc2_parameter("gamma",
                     function(n) runif(n, 0, 1),
                     function(x) dunif(x, 0, 1, log = TRUE),
                     min = 0, max = 1)))
  control <- smc2_control(n_parameter_sets, progress = FALSE)

  res <- smc2(p, f, control)
  expect_s3_class(res, "smc2_result")
  expect_setequal(names(res), c("pars", "probabilities", "statistics"))

  expect_equal(dim(res$pars), c(n_parameter_sets, 2))
  expect_equal(colnames(res$pars), c("beta", "gamma"))

  expect_equal(dim(res$probabilities), c(n_parameter_sets, 4))
  expect_equal(colnames(res$probabilities),
               c("log_prior", "log_likelihood", "log_posterior", "weight"))

  expect_setequal(
    names(res$statistics),
    c("ess", "acceptance_rate", "n_particles", "n_parameter_sets", "n_steps",
      "effort"))

  p <- predict.smc2_result(res)
  expect_equal(dim(p), dim(res$pars))
  expect_equal(dimnames(p), dimnames(res$pars))

  p <- predict.smc2_result(res, 100)
  expect_equal(dim(p), c(100, 2))
  expect_equal(dimnames(p), dimnames(res$pars))

  p <- predict.smc2_result(res, 1) # corner case
  expect_equal(dim(p), c(1, 2))
  expect_equal(dimnames(p), dimnames(res$pars))

  expect_gt(res$statistics$effort, 150)
  expect_lt(res$statistics$effort, 500)
})
