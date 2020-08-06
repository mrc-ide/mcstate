context("pmcmc")

## TODO: Ed and Lilith we could use better tests throughout here. The
## sampler should be ok to run for ~10k iterations without taking too
## long to be annoying in tests.
test_that("mcmc works for uniform distribution on unit square", {
  ## Uniform distribution:
  target <- function(p) 1

  proposal_kernel <- diag(2) * 0.1
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("a", "b")

  pars <- pmcmc_parameters$new(
    list(a = pmcmc_parameter(0.5, min = 0, max = 1),
         b = pmcmc_parameter(0.5, min = 0, max = 1)),
    proposal = proposal_kernel)

  set.seed(1)
  testthat::try_again(5, {
    res <- mcmc(pars, target, 1000, FALSE)
    expect_equal(res$acceptance_rate[["a"]], 1)
    expect_equal(res$acceptance_rate[["b"]], 1)
    expect_true(abs(mean(res$results$a) - 0.5) < 0.05)
    expect_true(abs(mean(res$results$b) - 0.5) < 0.05)
  })
})


test_that("mcmc works for multivariate gaussian", {
  target <- function(p) {
    mvtnorm::dmvnorm(unlist(p), log = TRUE)
  }

  proposal_kernel <- diag(2)
  pars <- pmcmc_parameters$new(
    list(a = pmcmc_parameter(0, min = -100, max = 100),
         b = pmcmc_parameter(0, min = -100, max = 100)),
    proposal = proposal_kernel)

  set.seed(1)
  testthat::try_again(5, {
    res <- mcmc(pars, target, 1000, FALSE)
    i <- seq(1, 1000, by = 20)
    expect_s3_class(res, "mcstate_pmcmc")
    expect_gt(ks.test(res$results$a[i], "pnorm")$p.value, 0.05)
    expect_gt(ks.test(res$results$b[i], "pnorm")$p.value, 0.05)
    expect_lt(abs(cov(res$results$a, res$results$b)), 0.1)
  })
})


test_that("multi-chain mcmc", {
  target <- function(p) 1

  proposal_kernel <- diag(2) * 0.1
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("a", "b")

  pars <- pmcmc_parameters$new(
    list(a = pmcmc_parameter(0.5, min = 0, max = 1),
         b = pmcmc_parameter(0.5, min = 0, max = 1)),
    proposal = proposal_kernel)

  res <- mcmc_multichain(pars, target, 1000, 3, FALSE)
})


test_that("pmcmc interface", {
  ## Mock particle filter:
  filter <- structure(list(run = target <- function(p) 1),
                      class = "particle_filter")

  proposal_kernel <- diag(2) * 0.1
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("a", "b")

  pars <- pmcmc_parameters$new(
    list(a = pmcmc_parameter(0.5, min = 0, max = 1),
         b = pmcmc_parameter(0.5, min = 0, max = 1)),
    proposal = proposal_kernel)

  set.seed(1)
  res <- pmcmc(pars, filter, 100)
  expect_is(res, "mcstate_pmcmc")
  set.seed(1)
  expect_identical(mcmc(pars, filter$run, 100, FALSE), res)
})


test_that("pmcmc passes expected arguments through to lower-level functions", {
  dat <- example_uniform()

  mock_mcmc <- mockery::mock()
  mock_mcmc_multichain <- mockery::mock()
  with_mock(
    "mcstate:::mcmc" = mock_mcmc,
    "mcstate:::mcmc_multichain" = mock_mcmc_multichain, {
      pmcmc(dat$pars, dat$filter, 100)
    })
  mockery::expect_called(mock_mcmc, 1)
  mockery::expect_called(mock_mcmc_multichain, 0)
  expect_equal(mockery::mock_args(mock_mcmc)[[1]],
               list(dat$pars, dat$target, 100, FALSE))
})


test_that("pmcmc can force multichain with just one chain", {
  dat <- example_uniform()

  mock_mcmc <- mockery::mock()
  mock_mcmc_multichain <- mockery::mock()
  with_mock(
    "mcstate:::mcmc" = mock_mcmc,
    "mcstate:::mcmc_multichain" = mock_mcmc_multichain, {
      pmcmc(dat$pars, dat$filter, 100, force_multichain = TRUE)
    })
  mockery::expect_called(mock_mcmc, 0)
  mockery::expect_called(mock_mcmc_multichain, 1)
  expect_equal(mockery::mock_args(mock_mcmc_multichain)[[1]],
               list(dat$pars, dat$target, 100, 1, FALSE))
})


test_that("can return proposals", {
  dat <- example_mvnorm()
  set.seed(1)
  res <- pmcmc(dat$pars, dat$filter, 50, return_proposals = TRUE)

  expect_s3_class(res$proposals, "data.frame")

  ## This is pretty lazy but shows that we're doing about the right thing
  p1 <- complex(real = res$results$a, imaginary = res$results$b)
  p2 <- complex(real = res$proposals$a, imaginary = res$proposals$b)
  expect_true(any(p2 %in% p1))
  expect_false(all(p2 %in% p1))

  ## Whenever we improved we definitely have the same answer:
  i <- res$results$log_likelihood[-1] > res$results$log_likelihood[-51]
  expect_equal(p2[which(i) + 1], p1[which(i) + 1])
})


test_that("can combine chains", {
  dat <- example_mvnorm()
  set.seed(1)
  res <- pmcmc(dat$pars, dat$filter, 100, return_proposals = TRUE,
               n_chains = 3)
  combined <- pmcmc_combine_chains(res, 51)

  expect_s3_class(combined, "data.frame")
  expect_equal(nrow(combined), 150)
  expect_equal(combined[1:50, ], rbind(res$chains[[1]]$results[52:101, ]),
               check.attributes = FALSE)
  expect_equal(combined[51:100, ], rbind(res$chains[[2]]$results[52:101, ]),
               check.attributes = FALSE)
  expect_equal(combined[101:150, ], rbind(res$chains[[3]]$results[52:101, ]),
               check.attributes = FALSE)
})


test_that("can't combine chains by taking more than burnin", {
  dat <- example_mvnorm()
  set.seed(1)
  res <- pmcmc(dat$pars, dat$filter, 100, n_chains = 3)
  expect_error(
    pmcmc_combine_chains(res, 101),
    "burn_in must be less than the total chain length")
  expect_equal(
    nrow(pmcmc_combine_chains(res, 100)),
    3)
})


test_that("Can't combine all chains with no burnin", {
  dat <- example_mvnorm()
  set.seed(1)
  res <- pmcmc(dat$pars, dat$filter, 10, n_chains = 3)
  expect_error(
    pmcmc_combine_chains(res, 0),
    "'burn_in' must be at least 1")
})


test_that("notify failure to compute gelman's diagnistic", {
  dat <- example_mvnorm()
  set.seed(1)
  expect_message(
    res <- pmcmc(dat$pars, dat$filter, 50, force_multichain = TRUE),
    "Could not calculate rhat: .+")
  expect_null(res$rhat)
})


test_that("Reflect parameters: double sided", {
  expect_equal(reflect_proposal(0.4, 0, 1), 0.4)
  expect_equal(reflect_proposal(1.4, 0, 1), 0.6)
  expect_equal(reflect_proposal(2.4, 0, 1), 0.4)

  ## Additional cases from the original
  expect_equal(reflect_proposal(6, 1, 5), 4)
  expect_equal(reflect_proposal(0, 1, 5), 2)
  expect_equal(reflect_proposal(10, 1, 5), 2)

  expect_equal(reflect_proposal(0:2 + 0.4, rep(0, 3), rep(1, 3)),
               c(0.4, 0.6, 0.4))
  expect_equal(reflect_proposal(0:2 + 1.4, rep(1, 3), rep(2, 3)),
               c(1.4, 1.6, 1.4))
})


test_that("Reflect parameters: infinite", {
  p <- rnorm(10)
  expect_equal(reflect_proposal(p, rep(-Inf, 5), rep(Inf, 5)), p)
})


test_that("reflect parameters: lower", {
  expect_equal(reflect_proposal(-0.4, 0, Inf), 0.4)
  expect_equal(reflect_proposal(0.6, 1, Inf), 1.4)
  expect_equal(reflect_proposal(0.6, 0, Inf), 0.6)
})


test_that("reflect parameters: upper", {
  expect_equal(reflect_proposal(1.4, -Inf, 1), 0.6)
  expect_equal(reflect_proposal(0.4, -Inf, 0), -0.4)
  expect_equal(reflect_proposal(0.4, -Inf, 1), 0.4)
})


test_that("run pmcmc with the particle filter", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  pars <- pmcmc_parameters$new(
    list(beta = pmcmc_parameter(0.2, min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         gamma = pmcmc_parameter(0.1, min = 0, max = 1,
                                 prior = function(p) log(1e-10))),
    proposal = proposal_kernel)

  dat <- example_sir()
  n_particles <- 100
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)

  ## A single chain
  set.seed(1)
  results <- pmcmc(pars, p, 30)
  expect_true(all(results$acceptance_rate[1:2] > 0))
})


test_that("proposal uses provided covariance structure", {
  ## Uniform distribution:
  target <- function(p) 1

  ## perfectly aligned (this is therefore not a valid mcmc for this
  ## distribution as it's not ergodic)
  proposal_kernel <- matrix(0.1, 2, 2)

  pars <- pmcmc_parameters$new(
    list(a = pmcmc_parameter(0.5, min = 0, max = 1),
         b = pmcmc_parameter(0.5, min = 0, max = 1)),
    proposal = proposal_kernel)

  res <- mcmc(pars, target, 100, FALSE)
  expect_equal(res$acceptance_rate[[1]], 1)
  expect_equal(res$results[[1]], res$results[[2]])
})
