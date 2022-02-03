test_that("Can run a tuning run", {
  dat <- example_mvnorm()
  pars <- dat$pars

  initial <- c(-5, -5)
  pars$update_proposal(diag(2) * 1e-3)

  ## Baseline control
  control <- pmcmc_control(100, save_trajectories = FALSE, save_state = FALSE,
                           progress = FALSE, n_chains = 3)
  res <- suppressMessages(
    pmcmc_tune(5, 100, pars = pars, filter = dat$filter, initial = initial,
               control = control))
  expect_setequal(names(res),
                  c("initial", "vcv", "proposal_kernel", "history"))
  expect_equal(res$proposal_kernel,
               r6_private(pars)$proposal_kernel)
  expect_equal(dim(res$initial), c(2, 3))
  expect_equal(dimnames(res$initial), list(c("a", "b"), NULL))
  expect_equal(res$vcv, res$proposal_kernel / (2.38^2 / 2))
  expect_setequal(names(res$history), c("pars", "probabilities"))
  expect_equal(dim(res$history$pars), c(2, 100, 3, 5))
  expect_equal(dim(res$history$probabilities), c(3, 100, 3, 5))
  expect_equal(rownames(res$history$pars), pars$names())
})


test_that("zero-length tuning does nothing", {
  dat <- example_mvnorm()
  pars <- dat$pars

  initial <- c(-5, -6)
  pars$update_proposal(diag(2) * 1e-3)

  ## Baseline control
  control <- pmcmc_control(100, save_trajectories = FALSE, save_state = FALSE,
                           progress = FALSE, n_chains = 3)
  res <- pmcmc_tune(0, 100, pars = pars, filter = dat$filter, initial = initial,
                    control = control)
  expect_equal(res$initial, c(a = -5, b = -6))
  expect_equal(res$vcv, matrix(NA_real_, 2, 2,
                               dimnames = list(c("a", "b"), c("a", "b"))))
  expect_equal(res$proposal_kernel, res$vcv)
  expect_equal(dim(res$history$pars), c(2, 100, 3, 0))
  expect_equal(dim(res$history$probabilities), c(3, 100, 3, 0))

  res2 <- pmcmc_tune(0, 100, pars = pars, filter = dat$filter,
                     control = control)
  expect_null(res2$initial)
  expect_equal(res2[names(res2) != "initial"], res[names(res) != "initial"])
})
