context("pMCMC")

test_that("MCMC can run", {
  range <- data.frame(name = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "model_data",
                      stringsAsFactors = FALSE)
  lprior <- list('beta' = function(pars) log(1e-10),
                 'gamma' = function(pars) log(1e-10))
  proposal_kernel <- diag(nrow(range)) * 0.01^2
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- range$name
  
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  n_particles <- 10
  n_mcmc <- 100
  
  # A single chain
  n_chains <- 1
  X <- pmcmc(range, lprior, p, n_particles, n_mcmc, proposal_kernel, 
             n_chains = n_chains)
  expect_equal(class(X), "mcstate_pmcmc")
  expect_equal(dim(X$results), c(n_mcmc + 1L, 5))
  expect_setequal(colnames(X$results), c(range$name, 'log_prior', 'log_likelihood', 'log_posterior'))
  expect_length(X$acceptance_rate, 5)
  expect_length(X$ess, 5)
  expect_true(all(X$acceptance_rate >= 0 & X$acceptance_rate <= 1))
  expect_true(all(X$ess >= 0))

  #plot(X)
  #summary(X)
  
  # Multiple chains
  n_chains <- 3
  burn_in <- 10
  Y <- pmcmc(range, lprior, p, n_particles, n_mcmc, proposal_kernel, 
             n_chains = n_chains)
  expect_equal(class(Y), "mcstate_pmcmc_list")
  expect_equal(nrow(Y$rhat$psrf), nrow(range))
  expect_length(Y$chains, 3)
  
  Y_master <- create_master_chain(Y, burn_in = burn_in)
  expect_equal(class(Y_master), "mcstate_pmcmc_chain")
  expect_length(Y_master, 5)
  expect_length(Y_master$beta, (n_mcmc + 1L - burn_in) * n_chains)

})

test_that("MCMC can infer the correct parameters", {
  skip("Test not yet complete")
  range <- data.frame(name = c("beta", "gamma"),
                      init = c(0.2, 0.1),
                      min = c(0, 0),
                      max = c(1, 1),
                      discrete = c(FALSE, FALSE),
                      target = "model_data",
                      stringsAsFactors = FALSE)
  lprior <- list('beta' = function(pars) log(1e-10),
                 'gamma' = function(pars) log(1e-10))
  proposal_kernel <- diag(nrow(range)) * 0.01^2
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- range$name
  
  dat <- example_sir()
  p <- particle_filter$new(dat$data, dat$model, dat$compare)
  n_particles <- 10
  n_mcmc <- 1000
  n_chains <- 1

  set.seed(1)
  X <- pmcmc(range, lprior, p, n_particles, n_mcmc, proposal_kernel, 
             n_chains = n_chains)

  expect_lt(abs(mean(X$results$beta[500:1000]) - 0.2), 0.01)
  expect_lt(abs(mean(X$results$gamma[500:1000]) - 0.1), 0.01)
})
