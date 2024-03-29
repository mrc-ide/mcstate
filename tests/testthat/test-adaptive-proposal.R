test_that("Can construct adaptive proposal", {
  control <- adaptive_proposal_control()
  pars <- example_sir()$pars
  obj <- adaptive_proposal$new(pars, control)
  expect_equal(obj$mean, pars$initial())
  expect_equal(obj$mean, pars$mean())
})


test_that("If mean is zero, autocorrelation is the vcv", {
  control <- adaptive_proposal_control()
  pars <- pmcmc_parameters$new(
    list(pmcmc_parameter("beta", 0.2, mean = 0),
         pmcmc_parameter("gamma", 0.1, mean = 0)),
    proposal = rbind(c(0.00057, 0.00034), c(0.00034, 0.00026)))

  obj <- adaptive_proposal$new(pars, control)
  expect_equal(unname(obj$autocorrelation), pars$vcv())
})

## initial autocor is zero if mean is zero

## weight needs to be at least 2 to start sensibly (w / (w - 1) is not
## well defined for w <= 1)

test_that("mean converges to weighted mean, regardless of acceptance", {
  control <- adaptive_proposal_control(initial_weight = 50)
  pars <- example_sir()$pars
  p <- pars$mean() * 2

  obj <- adaptive_proposal$new(pars, control)
  obj$proposal_was_adaptive <- TRUE
  for (i in 1:100) {
    obj$update(p, TRUE)
  }
  expect_equal(obj$mean, pars$mean() * 50 / 150 + p * 100 / 150)
  for (i in 1:100) {
    obj$update(p, FALSE)
  }
  expect_equal(obj$mean, pars$mean() * 50 / 250 + p * 200 / 250)
})


test_that("Scaling converges to expected limits", {
  control <- adaptive_proposal_control(initial_weight = 50)
  pars <- example_sir()$pars
  obj <- adaptive_proposal$new(pars, control)
  obj$proposal_was_adaptive <- TRUE
  p <- pars$mean()
  for (i in 1:100) {
    obj$update(p, TRUE)
  }
  expect_equal(
    obj$scaling,
    control$initial_scaling +
    100 * (1 - control$acceptance_target) * control$scaling_increment)

  obj <- adaptive_proposal$new(pars, control)
  obj$proposal_was_adaptive <- TRUE
  for (i in 1:100) {
    obj$update(p, FALSE)
  }
  expect_equal(
    obj$scaling,
    control$scaling_increment)
})
