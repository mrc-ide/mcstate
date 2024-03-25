test_that("Can construct adaptive proposal", {
  control <- adaptive_proposal_control()
  pars <- example_sir()$pars
  obj <- adaptive_proposal$new(pars, control)
  expect_equal(obj$mean, pars$initial())
  expect_equal(obj$mean, pars$mean())
})


test_that("mean and autocorrelation calculated correctly", {
  pars <- example_sir()$pars
  set.seed(1)
  n_pars <- 100
  d <- length(pars$mean())
  p <- t(replicate(n_pars, pars$propose(pars$mean())))
  p_list <- split(p, rep(seq_len(n_pars), d))
  
  control <- adaptive_proposal_control(initial_vcv_weight = 50,
                                       forget_rate = 0.1)
  obj <- adaptive_proposal$new(pars, control)
  for (i in 1:100) {
    obj$update(p[i, ], TRUE, p_list, i)
  }
  ## forget_rate = 0.1 so mean and autocorrelation should exclude first 10
  ## parameter sets
  expect_equal(obj$weight, 90)
  expect_equal(obj$included, seq(11, 100))
  expected_mean <- colMeans(p[11:100, ])
  expect_equal(obj$mean, expected_mean)
  expected_vcv <- cov(p[11:100, ])
  expect_equal(obj$vcv, expected_vcv)
  
  control <- adaptive_proposal_control(initial_vcv_weight = 50,
                                       forget_rate = 0.5,
                                       forget_end = 50)
  obj <- adaptive_proposal$new(pars, control)
  for (i in 1:100) {
    obj$update(p[i, ], TRUE, p_list, i)
  }
  ## forget_rate = 0.5 and forget_end = 50, so mean and autocorrelation should
  ## exclude first 25 parameter sets (half of the first 50)
  expect_equal(obj$weight, 75)
  expect_equal(obj$included, seq(26, 100))
  expected_mean <- colMeans(p[26:100, ])
  expect_equal(obj$mean, expected_mean)
  expected_vcv <- cov(p[26:100, ])
  expect_equal(obj$vcv, expected_vcv)
})


test_that("Scaling converges to expected limits - no diminishing adaptation", {
  control <- adaptive_proposal_control(initial_vcv_weight = 50,
                                       forget_rate = 0,
                                       pre_diminish = Inf)
  pars <- example_sir()$pars
  obj <- adaptive_proposal$new(pars, control)
  p <- pars$mean()
  for (i in 1:1000) {
    obj$update(p, 0.25, NULL, i)
  }
  expected_scaling <-
    exp(log(control$initial_scaling) + 1000 * (0.25 - control$acceptance_target)
        * obj$scaling_increment / sqrt(obj$scaling_weight))
  expect_equal(obj$scaling, expected_scaling) 

  control <- adaptive_proposal_control(initial_vcv_weight = 50,
                                       forget_rate = 0,
                                       pre_diminish = Inf,
                                       log_scaling_update = FALSE)
  pars <- example_sir()$pars
  obj <- adaptive_proposal$new(pars, control)
  p <- pars$mean()
  for (i in 1:1000) {
    obj$update(p, 0.25, NULL, i)
  }
  expected_scaling <-
    control$initial_scaling + 1000 * (0.25 - control$acceptance_target) * 
    obj$scaling_increment / sqrt(obj$scaling_weight)
  expect_equal(obj$scaling, expected_scaling) 
})


test_that("Scaling converges to expected limits - diminishing adaptation", {
  control <- adaptive_proposal_control(initial_vcv_weight = 50,
                                       forget_rate = 0,
                                       pre_diminish = 0)
  pars <- example_sir()$pars
  obj <- adaptive_proposal$new(pars, control)
  initial_scaling_weight <- obj$scaling_weight
  p <- pars$mean()
  for (i in 1:1000) {
    obj$update(p, 0.25, NULL, i)
  }
  expected_scaling <-
    exp(log(control$initial_scaling) + (0.25 - control$acceptance_target)
        * sum(obj$scaling_increment / 
                sqrt(initial_scaling_weight + seq_len(1000))))
  expect_equal(obj$scaling, expected_scaling) 
  
  control <- adaptive_proposal_control(initial_vcv_weight = 50,
                                       forget_rate = 0,
                                       pre_diminish = 0,
                                       log_scaling_update = FALSE)
  pars <- example_sir()$pars
  obj <- adaptive_proposal$new(pars, control)
  initial_scaling_weight <- obj$scaling_weight
  p <- pars$mean()
  for (i in 1:1000) {
    obj$update(p, 0.25, NULL, i)
  }
  expected_scaling <-
    control$initial_scaling + (0.25 - control$acceptance_target) * 
    sum(obj$scaling_increment / sqrt(initial_scaling_weight + seq_len(1000)))
  expect_equal(obj$scaling, expected_scaling)
})
