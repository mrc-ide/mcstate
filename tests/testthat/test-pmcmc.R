context("pmcmc")

test_that("MCMC can run", {
  skip("old test")
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  pars <- pmcmc_parameters(
    beta = pmcmc_parameter(0.2, min = 0, max = 1,
                           prior = function(p) log(1e-10)),
    gamma = pmcmc_parameter(0.1, min = 0, max = 1,
                            prior = function(p) log(1e-10)),
    .proposal = proposal_kernel)

  dat <- example_sir()
  n_particles <- 10
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)

  n_mcmc <- 100

  # A single chain
  n_chains <- 1
  mcmc_results <- pmcmc(pars, p, n_mcmc, n_chains = n_chains)

  expect_equal(class(mcmc_results), "mcstate_pmcmc")
  expect_equal(dim(mcmc_results$results), c(n_mcmc + 1L, 5))
  expect_setequal(colnames(mcmc_results$results),
                  c(range$name,
                    "log_prior", "log_likelihood", "log_posterior"))
  expect_length(mcmc_results$acceptance_rate, 5)
  expect_length(mcmc_results$ess, 5)
  expect_true(all(mcmc_results$acceptance_rate >= 0 &
                  mcmc_results$acceptance_rate <= 1))
  expect_true(all(mcmc_results$ess >= 0))

  summary(mcmc_results)
  plot(mcmc_results)

  # Multiple chains
  n_chains <- 3
  burn_in <- 10
  multi_chain <- pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                       n_chains = n_chains)

  expect_equal(class(multi_chain), "mcstate_pmcmc_list")
  expect_equal(nrow(multi_chain$rhat$psrf), nrow(range))
  expect_length(multi_chain$chains, 3)

  summary(multi_chain)
  plot(multi_chain)

  multi_chain_master <- create_master_chain(multi_chain, burn_in = burn_in)
  expect_equal(ncol(multi_chain_master), 5)
  expect_equal(nrow(multi_chain_master), (n_mcmc + 1L - burn_in) * n_chains)

  summary(multi_chain_master)
  plot(multi_chain_master)
})

test_that("MCMC doesn't move away from correct parameters", {
  skip("old test")
  range <- data.frame(name = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  lprior <- list("beta" = function(pars) log(1e-10),
                 "gamma" = function(pars) log(1e-10))
  proposal_kernel <- diag(nrow(range)) * 0.01^2
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- range$name

  dat <- example_sir()
  n_particles <- 20
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  n_mcmc <- 500
  n_chains <- 1

  set.seed(1)
  mcmc_results <- pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                        n_chains = n_chains)

  expect_lt(abs(mean(mcmc_results$results$beta) - 0.2), 0.03)
  expect_lt(abs(mean(mcmc_results$results$gamma) - 0.1), 0.02)
})

test_that("MCMC runs on different targets", {
  skip("old test")
  # All three targets
  range <- data.frame(name = c("beta", "step_start", "exp_noise"),
                      init = c(0.2, 200, 1e6),
                      min = c(0, 0, 1),
                      max = c(1, 400, 1e10),
                      discrete = c(FALSE, TRUE, FALSE),
                      target = c("pars_model", "pars_initial", "pars_compare"),
                      stringsAsFactors = FALSE)
  # One informative prior
  lprior <- list("beta" = function(pars) dnorm(pars["beta"], 0.2, 0.01),
                 "step_start" = function(pars) log(1e-10),
                 "exp_noise" = function(pars) log(1e-10))
  proposal_kernel <- matrix(c(0.001^2, 0, 0,
                              0, 0.001^2, 0,
                              0,       0, 10^2),
                            nrow = 3, byrow = TRUE,
                            dimnames = list(
                              range$name,
                              range$name))
  initial <- function(info, n_particles, pars) {
    list(step = pars[["step_start"]])
  }

  dat <- example_sir()
  data <- dat$data
  offset <- 400
  data[c("step_start", "step_end")] <-
    data[c("step_start", "step_end")] + offset
  data$step_start[[1]] <- 0
  n_particles <- 10
  p <- particle_filter$new(data, dat$model, n_particles, dat$compare,
                           index = dat$index, initial = initial)
  n_mcmc <- 10

  # A single chain
  n_chains <- 1
  mcmc_results <- pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                        n_chains = n_chains)

  expect_equal(class(mcmc_results), "mcstate_pmcmc")
  expect_equal(dim(mcmc_results$results), c(n_mcmc + 1L, 6))
  expect_setequal(colnames(mcmc_results$results),
                  c(range$name,
                    "log_prior", "log_likelihood", "log_posterior"))
  expect_length(mcmc_results$acceptance_rate, 6)
  expect_length(mcmc_results$ess, 6)
  expect_true(all(mcmc_results$acceptance_rate >= 0 &
                  mcmc_results$acceptance_rate <= 1))
  expect_true(all(mcmc_results$ess >= 0))

  summary(mcmc_results)
  plot(mcmc_results)
})

test_that("Proposals are correctly reflected", {
  skip("old test")
  expect_equal(object = reflect_proposal(x = 6, floor = 1, cap = 5),
               expected = 4)
  expect_equal(object = reflect_proposal(x = 0, floor = 1, cap = 5),
               expected = 2)
  expect_equal(object = reflect_proposal(x = 10, floor = 1, cap = 5),
               expected = 2)

  # check that the function behaves as expected when passed a vector
  n <- 10
  tmp <- data.frame(x = rnorm(n, 1),
                    floor = runif(n))
  tmp$cap <- tmp$floor + runif(n)

  df <- with(tmp, reflect_proposal(x, floor, cap))
  vec <- with(tmp, mapply(FUN = reflect_proposal,
                          x = x,
                          floor = floor,
                          cap = cap))

  expect_equal(df, vec)
})

test_that("MCMC jumps behave as expected", {
  skip("old test")
  ## check that proposing jumps of size zero results in the
  ## initial parameter being retained
  range <- data.frame(name = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  lprior <- list("beta" = function(pars) log(1e-10),
                 "gamma" = function(pars) log(1e-10))
  proposal_kernel <- matrix(rep(0, 4),
                            nrow = 2, byrow = TRUE,
                            dimnames = list(
                              range$name,
                              range$name))

  dat <- example_sir()
  n_particles <- 20
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  n_mcmc <- 500
  n_chains <- 1
  mcmc_results <- pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                        n_chains = n_chains, output_proposals = TRUE)

  expect_equal(object = mcmc_results$results$beta,
               expected = rep(range$init[1], n_mcmc + 1))
  expect_equal(object = mcmc_results$results$gamma,
               expected = rep(range$init[2], n_mcmc + 1))
  expect_true(!all(diff(mcmc_results$results$log_likelihood) == 0))
  expect_equal(dim(mcmc_results$proposals), c(n_mcmc + 1L, 6))

  ## check non-zero covariance ihas an impact on proposals
  set.seed(1)
  proposal_kernel <- matrix(c(0.01^2, 0,
                              0, 0.01^2),
                            nrow = 2, byrow = TRUE,
                            dimnames = list(
                              range$name,
                              range$name))
  no_covar_results <- pmcmc(range, lprior, p, n_mcmc,
                            proposal_kernel, n_chains = n_chains,
                            output_proposals = TRUE)
  set.seed(1)
  proposal_kernel <- matrix(c(0.01^2, 0.01^2,
                              0.01^2, 0.01^2),
                            nrow = 2, byrow = TRUE,
                            dimnames = list(
                              range$name,
                              range$name))
  covar_results <- pmcmc(range, lprior, p, n_mcmc,
                         proposal_kernel, n_chains = n_chains,
                         output_proposals = TRUE)
  expect_false(all(covar_results$results == no_covar_results$results))
})

test_that("MCMC range input errors", {
  skip("old test")
  range <- data.frame(name = c("beta", "gamma"),
                      init = c(-0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  expect_error(mcmc_validate_range(range),
               "initial parameters are outside of specified range")
  range <- data.frame(name = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(2, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  expect_error(mcmc_validate_range(range),
               "'discrete' entries must be TRUE or FALSE")
  range <- data.frame(name = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c("0", 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  expect_error(mcmc_validate_range(range),
               "'min' entries must be numeric")
  range <- data.frame(name = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c("1", 1),
                      discrete = c(FALSE, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  expect_error(mcmc_validate_range(range),
               "'max' entries must be numeric")
  range <- data.frame(name = c("beta", "gamma"),
                      init = c("0.2", 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  expect_error(mcmc_validate_range(range),
               "'init' entries must be numeric")
  range <- data.frame(name = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "model_pars",
                      stringsAsFactors = FALSE)
  expect_error(mcmc_validate_range(range),
               paste0("Invalid target 'model_pars': must be one of ",
                      "'pars_model', 'pars_compare', 'pars_initial'"))
  range <- data.frame(names = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  expect_error(mcmc_validate_range(range),
               "Missing columns from range: 'name'")
  range <- data.frame(name = c("beta", "beta"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  expect_error(mcmc_validate_range(range),
               "Duplicate 'name' entries not allowed in range")
})

test_that("MCMC function input errors", {
  skip("old test")
  range <- data.frame(name = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  lprior <- list("beta" = function(pars) log(1e-10),
                 "gamma" = function(pars) log(1e-10))
  proposal_kernel <- diag(nrow(range)) * 0.01^2
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- range$name

  dat <- example_sir()
  n_particles <- 10
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  n_mcmc <- 100

  n_chains <- 1
  expect_error(pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                     n_chains = n_chains, output_proposals = "yes_please"),
               "output_proposals must be either TRUE or FALSE")

  # bad priors
  lprior <- list("beta" = function(pars) log(1e-10))
  expect_error(pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                     n_chains = n_chains),
               "All sampled parameters must have a defined prior")
  lprior <- list("beta" = function(pars) c(log(1e-10), log(1e-10)),
                 "gamma" = function(pars) c(0, 0))
  expect_error(pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                     n_chains = n_chains),
               paste0("lprior_funcs must return a single numeric representing ",
                      "the log prior"))
  lprior <- list("beta" = function(pars) log(0),
                 "gamma" = function(pars) 0)
  expect_error(pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                     n_chains = n_chains),
               "initial parameters are not compatible with supplied prior")

  # incomplete proposal kernel
  lprior <- list("beta" = function(pars) log(1e-10),
                 "gamma" = function(pars) log(1e-10))
  proposal_kernel <- diag(1) * 0.01^2
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- range$name[1]
  expect_error(pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                     n_chains = n_chains),
               paste0("proposal_kernel must be a matrix or vector with names ",
                      "corresponding to the parameters being sampled"))
})

test_that("Master chain errors", {
  skip("old test")
  expect_error(create_master_chain(list(rep(NA, 10))),
               "x must be a pmcmc_list object")


  range <- data.frame(name = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "pars_model",
                      stringsAsFactors = FALSE)
  lprior <- list("beta" = function(pars) log(1e-10),
                 "gamma" = function(pars) log(1e-10))
  proposal_kernel <- diag(nrow(range)) * 0.01^2
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- range$name

  dat <- example_sir()
  n_particles <- 10
  p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                           index = dat$index)
  n_mcmc <- 10
  n_chains <- 3
  multi_chain <- pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                       n_chains = n_chains)

  expect_error(create_master_chain(multi_chain, burn_in = -1L),
               "burn_in must not be negative")
  expect_error(create_master_chain(multi_chain, burn_in = "10"),
               "burn_in must be an integer")
  expect_error(create_master_chain(multi_chain, burn_in = 11),
               "burn_in is greater than chain length")
})


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


test_that("can compute summary of a chain", {
  dat <- example_mvnorm()
  set.seed(1)
  res <- pmcmc(dat$pars, dat$filter, 50, return_proposals = TRUE)
  ans <- summary(res)
  expect_type(ans, "list")
  expect_setequal(names(ans), c("summary", "corr_mat", "sd"))
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
