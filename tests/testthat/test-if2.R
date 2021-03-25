context("if2")


test_that("Can run IF2", {
  pars <- list(beta = 0.15, gamma = 0.05)
  pars_sd <- c(beta = 0.02, gamma = 0.02) # not a list!

  dat <- example_sir()

  iterations <- 50
  cooling_target <- 0.5
  n_par_sets <- 100
  set.seed(1)
  if2_res <- if2(pars, dat$data, dat$model, dat$compare, NULL, dat$index,
                 pars_sd, iterations, n_par_sets, cooling_target)

  expect_setequal(
    names(if2_res),
    c("log_likelihood", "if_pars"))

  expect_equal(dim(if2_res$log_likelihood), c(iterations, n_par_sets))
  expect_equal(dim(if2_res$if_pars),
               c(length(pars), iterations, n_par_sets))
})

