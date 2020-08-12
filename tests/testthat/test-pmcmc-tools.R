context("pmcmc (tools)")

test_that("pmcmc_thin with no args is a no-op", {
  results <- example_sir_pmcmc()$pmcmc
  expect_identical(pmcmc_thin(results), results)
})


test_that("discarding burnin drops beginnings of chain", {
  results <- example_sir_pmcmc()$pmcmc
  res <- pmcmc_thin(results, 10)
  i <- 11:31
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
  expect_identical(res$state, results$state[, i])
  expect_identical(res$trajectories, results$trajectories[, , i])
})


test_that("thinning drops all over chain", {
  results <- example_sir_pmcmc()$pmcmc
  res <- pmcmc_thin(results, thin = 4)
  i <- seq(1, 31, by = 4)
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
  expect_identical(res$state, results$state[, i])
  expect_identical(res$trajectories, results$trajectories[, , i])
})


test_that("burnin and thin can be used together", {
  results <- example_sir_pmcmc()$pmcmc
  i <- seq(11, 31, by = 4)
  res <- pmcmc_thin(results, 10, 4)
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
  expect_identical(res$state, results$state[, i])
  expect_identical(res$trajectories, results$trajectories[, , i])
})


test_that("can discard whole chain", {
  results <- example_sir_pmcmc()$pmcmc
  res <- pmcmc_thin(results, 31)
  i <- integer(0)
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
  expect_identical(res$state, results$state[, i])
  expect_identical(res$trajectories, results$trajectories[, , i])
})


test_that("can't discard more than the whole chain", {
  results <- example_sir_pmcmc()$pmcmc
  expect_error(pmcmc_thin(results, 32),
               "'burnin' can be at most 31 for your results")
})


test_that("Can thin when no state/trajectories present", {
  results <- example_sir_pmcmc()$pmcmc
  results$trajectories <- NULL
  results$state <- NULL

  i <- seq(11, 31, by = 4)
  res <- pmcmc_thin(results, 10, 4)
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
  expect_null(res$state)
  expect_null(res$trajectories)
})
