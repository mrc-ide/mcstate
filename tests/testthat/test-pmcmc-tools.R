context("pmcmc (tools)")

test_that("pmcmc_thin with no args is a no-op", {
  results <- example_sir_pmcmc()$pmcmc
  expect_identical(pmcmc_thin(results), results)
})


test_that("discarding burnin drops beginnings of chain", {
  results <- example_sir_pmcmc()$pmcmc
  res <- pmcmc_thin(results, 9)
  i <- 10:30
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
  expect_identical(res$state, results$state[, i])
  expect_identical(res$trajectories$state, results$trajectories$state[, i, ])
  expect_identical(res$restart$state,
                   results$restart$state[, i, , drop = FALSE])
})


test_that("thinning drops all over chain", {
  results <- example_sir_pmcmc()$pmcmc
  res <- pmcmc_thin(results, thin = 4)
  i <- seq(1, 30, by = 4)
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
  expect_identical(res$state, results$state[, i])
  expect_identical(res$trajectories$state, results$trajectories$state[, i, ])
  expect_identical(res$restart$state,
                   results$restart$state[, i, , drop = FALSE])
})


test_that("burnin and thin can be used together", {
  results <- example_sir_pmcmc()$pmcmc
  i <- seq(10, 30, by = 4)
  res <- pmcmc_thin(results, 9, 4)
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
  expect_identical(res$state, results$state[, i])
  expect_identical(res$trajectories$state, results$trajectories$state[, i, ])
  expect_identical(res$restart$state,
                   results$restart$state[, i, , drop = FALSE])
})


test_that("can't discard the whole chain (or more)", {
  results <- example_sir_pmcmc()$pmcmc
  expect_error(pmcmc_thin(results, 30),
               "'burnin' must be less than 30 for your results")
  expect_error(pmcmc_thin(results, 100),
               "'burnin' must be less than 30 for your results")
})


test_that("Can thin when no state/trajectories present", {
  results <- example_sir_pmcmc()$pmcmc
  results$trajectories <- NULL
  results$state <- NULL

  i <- seq(10, 30, by = 4)
  res <- pmcmc_thin(results, 9, 4)
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
  expect_null(res$state)
  expect_null(res$trajectories)
})


test_that("can combine chains", {
  results <- example_sir_pmcmc2()$results

  results1 <- results[[1]]
  results2 <- results[[2]]
  results3 <- results[[3]]

  res <- pmcmc_combine(results1, results2, results3)

  n_steps <- nrow(results1$pars)
  n_par <- ncol(results1$pars)
  n_particles <- nrow(results1$state)
  n_index <- nrow(results1$trajectories$state)
  n_time <- dim(results1$trajectories$state)[[3]]
  n_restart <- dim(results1$restart$state)[[3]]
  n_state <- nrow(results1$state)

  n_steps3 <- n_steps * 3

  expect_equal(dim(res$pars), c(n_steps3, n_par))
  expect_equal(dim(res$probabilities), c(n_steps3, 3))
  expect_equal(dim(res$state), c(n_state, n_steps3))
  expect_equal(dim(res$trajectories$state), c(n_index, n_steps3, n_time))
  expect_equal(dim(res$restart$state), c(n_state, n_steps3, n_restart))

  i <- seq_len(n_steps) + n_steps
  expect_equal(res$pars[i, ], results2$pars)
  expect_equal(res$probabilities[i, ], results2$probabilities)
  expect_equal(res$state[, i], results2$state)
  expect_equal(res$trajectories$state[, i, ], results2$trajectories$state)
  expect_equal(res$restart$state[, i, , drop = FALSE], results2$restart$state)
})


test_that("can combine chains with list interface", {
  results <- example_sir_pmcmc2()$results
  expect_identical(
    pmcmc_combine(results[[1]], results[[2]], results[[3]]),
    pmcmc_combine(samples = results))
})


test_that("can combine chains without samples or state", {
  f <- function(x) {
    x$state <- NULL
    x$trajectories <- NULL
    x$predict <- NULL
    x
  }

  results1 <- example_sir_pmcmc2()$results
  results2 <- lapply(results1, f)
  combined1 <- pmcmc_combine(samples = results1)
  combined2 <- pmcmc_combine(samples = results2)

  expect_null(combined2$state)
  expect_null(combined2$trajectories)
  expect_null(combined2$predict)

  nms <- setdiff(names(combined1), c("state", "trajectories", "predict"))
  expect_identical(combined1[nms], combined2[nms])
})


test_that("can drop burnin from combined chains", {
  results <- example_sir_pmcmc2()$results
  combined <- pmcmc_combine(samples = results)
  res <- pmcmc_thin(combined, burnin = 9)
  expect_equal(res$chain, rep(1:3, each = 21))
  expect_equal(res$iteration, rep(10:30, 3))

  ## Same performed either way:
  expect_identical(
    res,
    pmcmc_combine(samples = lapply(results, pmcmc_thin, burnin = 9)))
})


test_that("can thin combined chains", {
  results <- example_sir_pmcmc2()$results
  combined <- pmcmc_combine(samples = results)
  res <- pmcmc_thin(combined, burnin = 9, thin = 4)
  expect_equal(res$chain, rep(1:3, each = 6))
  expect_equal(res$iteration, rep(seq(10, 30, by = 4), 3))

  ## Same performed either way:
  expect_identical(
    res,
    pmcmc_combine(samples = lapply(results, pmcmc_thin, 9, 4)))
})


test_that("combining ignores names", {
  results <- example_sir_pmcmc2()$results
  expect_identical(
    pmcmc_combine(samples = results),
    pmcmc_combine(samples = list(a = results[[1]],
                                 b = results[[2]],
                                 c = results[[3]])))
  expect_identical(
    pmcmc_combine(samples = results),
    pmcmc_combine(a = results[[1]], b = results[[2]], c = results[[3]]))
})


test_that("combining requires at least one run", {
  results <- example_sir_pmcmc2()$results
  expect_error(pmcmc_combine(),
               "At least 2 samples objects must be provided")
  expect_error(pmcmc_combine(samples = NULL),
               "'samples' must be a list")
  expect_error(pmcmc_combine(results[[1]]),
               "At least 2 samples objects must be provided")
})


test_that("can't recombine chains", {
  results <- example_sir_pmcmc2()$results
  expect_error(
    pmcmc_combine(results[[1]], pmcmc_combine(results[[2]], results[[3]])),
    "Chains have already been combined")
})


test_that("require consistent data", {
  results <- example_sir_pmcmc2()$results
  a <- results[[1]]
  b <- results[[2]]
  expect_error(
    pmcmc_combine(a, pmcmc_thin(b, burnin = 2)),
    "All chains must have the same length")
})


test_that("can't combine chains with different parameters", {
  results <- example_sir_pmcmc2()$results
  a <- results[[1]]
  a$pars <- cbind(a$pars, zeta = 1)
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "All parameters must have the same dimension names")
})


test_that("combine require the same iterations", {
  results <- example_sir_pmcmc2()$results
  a <- results[[1]]
  a$iteration <- a$iteration + 1
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "All chains must have the same iterations")
})


test_that("Can't combine chains that differ in if they have state", {
  results <- example_sir_pmcmc2()$results
  a <- results[[1]]
  a$state <- NULL
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "If 'state' is present for any samples, it must be present for all")
})


test_that("Can't combine chains that differ in if they have trajectories", {
  results <- example_sir_pmcmc2()$results
  a <- results[[1]]
  a$trajectories <- NULL
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "If 'trajectories' is present for any samples, it must be present for all")
})


test_that("Can't combine chains that differ in if they have restart", {
  results <- example_sir_pmcmc2()$results
  a <- results[[1]]
  a$restart <- NULL
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "If 'restart' is present for any samples, it must be present for all")
})


test_that("Can't combine inconsistent trajectories", {
  results <- example_sir_pmcmc2()$results
  a <- results[[1]]
  a$trajectories$rate <- a$trajectories$rate + 1
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "trajectories data is inconsistent")
})


test_that("check object types for combine", {
  results <- example_sir_pmcmc2()$results
  expect_error(
    pmcmc_combine(results[[1]], NULL),
    "Elements of", fixed = TRUE)
  expect_error(
    pmcmc_combine(samples = list(results[[1]], NULL)),
    "Elements of", fixed = TRUE)
})


test_that("can sample from a mcmc", {
  results <- example_sir_pmcmc()$pmcmc
  sub <- pmcmc_sample(results, 10, burnin = 9)
  expect_equal(nrow(sub$pars), 10)
  expect_true(all(sub$iteration >= 10))
})


test_that("sampling is with replacement", {
  results <- example_sir_pmcmc()$pmcmc
  sub <- pmcmc_sample(results, 50, burnin = 9)
  expect_equal(nrow(sub$pars), 50)
  expect_true(all(sub$iteration >= 10))
  expect_true(any(duplicated(sub$iteration)))
})


test_that("can sample from a combined chain", {
  results <- pmcmc_combine(samples = example_sir_pmcmc2()$results)
  sub <- pmcmc_sample(results, 50, burnin = 9)
  expect_equal(nrow(sub$pars), 50)
  expect_true(all(1:3 %in% sub$chain))
  expect_true(all(sub$iteration >= 10))
})


test_that("combining chains keeps rownames", {
  results <- example_sir_pmcmc2()$results

  nms_t <- c("S", "I", "R")
  nms_s <- c("S", "I", "R", "cum", "inc")
  for (i in seq_along(results)) {
    rownames(results[[i]]$trajectories$state) <- nms_t
    rownames(results[[i]]$state) <- nms_s
  }

  results1 <- results[[1]]
  results2 <- results[[2]]
  results3 <- results[[3]]

  res <- pmcmc_combine(results1, results2, results3)

  expect_equal(rownames(res$state), nms_s)
  expect_equal(rownames(res$trajectories$state), nms_t)
})

test_that("Can't combine inconsistent nested trajectories", {
  results <- example_sir_nested_pmcmc()$results
  a <- results[[1]]
  a$trajectories$rate <- a$trajectories$rate + 1
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "trajectories data is inconsistent")
})

test_that("can't recombine nested chains", {
  results <- example_sir_nested_pmcmc()$results
  expect_error(
    pmcmc_combine(results[[1]],
                         pmcmc_combine(results[[2]], results[[3]])),
    "Chains have already been combined")
})

test_that("nested combining requires at least one run", {
  results <- example_sir_nested_pmcmc()$results
  expect_error(pmcmc_combine(),
               "At least 2 samples objects must be provided")
  expect_error(pmcmc_combine(samples = list()),
               "At least 2 samples objects must be provided")
  expect_error(pmcmc_combine(results[[1]]),
               "At least 2 samples objects must be provided")
})

test_that("can't combine nested chains with different parameters", {
  results <- example_sir_nested_pmcmc()$results
  a <- results[[1]]
  a$pars <- cbind(a$pars, zeta = 1)
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "All parameters must have the same dimension names")
})


test_that("nested example_sir_nested_pmcmc require the same iterations", {
  results <- example_sir_nested_pmcmc()$results
  a <- results[[1]]
  a$iteration <- a$iteration + 1
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "All chains must have the same iterations")
})

test_that("require consistent nested data", {
  results <- example_sir_nested_pmcmc()$results
  a <- results[[1]]
  b <- results[[2]]
  expect_error(
    pmcmc_combine(a, pmcmc_thin(b, burnin = 2)),
    "All chains must have the same length")
})


test_that("discarding burnin drops beginnings of nested chain", {
  results <- example_sir_nested_pmcmc()$results[[1]]
  res <- pmcmc_thin(results, 9)
  i <- 10:30
  expect_identical(res$pars, results$pars[i, , ])
  expect_identical(res$probabilities, results$probabilities[i, , ])
  expect_identical(res$state, results$state[, , i])
  expect_identical(res$trajectories$state, results$trajectories$state[, , i, ])
  expect_identical(res$restart$state,
                   results$restart$state[, , i, , drop = FALSE])
})

test_that("can sample from a nested mcmc", {
  results <- example_sir_nested_pmcmc()$results[[1]]
  sub <- pmcmc_sample(results, 10, burnin = 9)
  expect_equal(nrow(sub$pars), 10)
  expect_true(all(sub$iteration >= 10))
})

test_that("Can't combine nested chains that differ in if they have state", {
  results <- example_sir_nested_pmcmc()$results
  a <- results[[1]]
  a$state <- NULL
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "If 'state' is present for any samples, it must be present for all")
})

test_that("Can combine nested chains without save state", {
  results <- example_sir_nested_pmcmc()$results
  results[[1]]$state <- NULL
  results[[2]]$state <- NULL
  expect_equal(pmcmc_combine(results[[1]], results[[2]])$state, NULL)
})

test_that("Can't combine nested chains that differ if have trajectories", {
  results <- example_sir_nested_pmcmc()$results
  a <- results[[1]]
  a$trajectories <- NULL
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "If 'trajectories' is present for any samples, it must be present for all")
})


test_that("Can't combine nested chains that differ in if they have restart", {
  results <- example_sir_nested_pmcmc()$results
  a <- results[[1]]
  a$restart <- NULL
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "If 'restart' is present for any samples, it must be present for all")
})


test_that("Can't combine inconsistent nested trajectories", {
  results <- example_sir_nested_pmcmc()$results
  a <- results[[1]]
  a$trajectories$rate <- a$trajectories$rate + 1
  expect_error(
    pmcmc_combine(a, results[[2]]),
    "trajectories data is inconsistent")
})


test_that("can combine chains for nested model", {
  results <- example_sir_nested_pmcmc()$results

  results1 <- results[[1]]
  results2 <- results[[2]]
  results3 <- results[[3]]

  res <- pmcmc_combine(results1, results2, results3)

  n_steps <- nrow(results1$pars)
  n_par <- ncol(results1$pars)
  n_pop <- nlayer(results1$pars)
  n_particles <- nrow(results1$state)
  n_index <- nrow(results1$trajectories$state)
  n_time <- dim(results1$trajectories$state)[[4]]
  n_restart <- dim(results1$restart$state)[[4]]
  n_state <- nrow(results1$state)

  n_steps3 <- n_steps * 3

  expect_equal(dim(res$pars), c(n_steps3, n_par, n_pop))
  expect_equal(dim(res$probabilities), c(n_steps3, 3, n_pop))
  expect_equal(dim(res$state), c(n_state, n_pop, n_steps3))
  expect_equal(dim(res$trajectories$state), c(n_index, n_pop, n_steps3, n_time))
  expect_equal(dim(res$restart$state), c(n_state, n_pop, n_steps3, n_restart))

  i <- seq_len(n_steps) + n_steps
  expect_equal(res$pars[i, , ], results2$pars)
  expect_equal(res$probabilities[i, , ], results2$probabilities)
  expect_equal(res$state[, , i], results2$state)
  expect_equal(res$trajectories$state[, , i, ], results2$trajectories$state)
  expect_equal(res$restart$state[, , i, , drop = FALSE],
               results2$restart$state)
})
