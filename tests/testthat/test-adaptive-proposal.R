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
  vcv <- obj$autocorrelation - obj$weight / (obj$weight - 1) * qp(obj$mean)
  expect_equal(vcv, expected_vcv)
  expect_equal(adaptive_vcv(0.5, obj$autocorrelation, obj$weight, obj$mean,
                            pars$vcv(), obj$control$initial_vcv_weight),
               0.5 * (89 * vcv + (50 + d + 1) * pars$vcv()) / (90 + 50 + d + 1))
  
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
  vcv <- obj$autocorrelation - obj$weight / (obj$weight - 1) * qp(obj$mean)
  expect_equal(vcv, expected_vcv)
  expect_equal(adaptive_vcv(0.5, obj$autocorrelation, obj$weight, obj$mean,
                            pars$vcv(), obj$control$initial_vcv_weight),
               0.5 * (74 * vcv + (50 + d + 1) * pars$vcv()) / (75 + 50 + d + 1))
})


test_that("Scaling converges to expected limits - no diminishing adaptation", {
  control <- adaptive_proposal_control(initial_vcv_weight = 50,
                                       forget_rate = 0,
                                       pre_diminish = Inf)
  pars <- example_sir()$pars
  obj <- adaptive_proposal$new(pars, control)
  p <- pars$mean()
  for (i in 1:1000) {
    obj$update(p, TRUE, NULL, i)
  }
  expect_equal(
    obj$scaling,
    (sqrt(control$initial_scaling) +
     1000 * (1 - control$acceptance_target) * control$scaling_increment) ^ 2) 

  obj <- adaptive_proposal$new(pars, control)
  for (i in 1:1000) {
    obj$update(p, FALSE, NULL, i)
  }
  expect_equal(
    obj$scaling,
    control$scaling_increment ^ 2)
})


test_that("Scaling converges to expected limits - diminishing adaptation", {
  control <- adaptive_proposal_control(initial_vcv_weight = 50,
                                       forget_rate = 0,
                                       pre_diminish = 0)
  pars <- example_sir()$pars
  obj <- adaptive_proposal$new(pars, control)
  p <- pars$mean()
  for (i in 1:1000) {
    obj$update(p, TRUE, NULL, i)
  }
  expect_equal(
    obj$scaling,
    (sqrt(control$initial_scaling) + sum(1 / sqrt(seq_len(1000))) *
       (1 - control$acceptance_target) * control$scaling_increment) ^ 2) 
  
  obj <- adaptive_proposal$new(pars, control)
  for (i in 1:1000) {
    obj$update(p, FALSE, NULL, i)
  }
  expect_equal(
    obj$scaling,
    (sqrt(control$initial_scaling) - sum(1 / sqrt(seq_len(1000))) *
      control$acceptance_target * control$scaling_increment) ^ 2)
})
