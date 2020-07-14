context("pmcmc")

test_that("MCMC can run", {
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

  # A single chain
  n_chains <- 1
  mcmc_results <- pmcmc(range, lprior, p, n_mcmc, proposal_kernel,
                        n_chains = n_chains)

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
