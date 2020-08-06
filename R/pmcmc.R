##' Run a pmcmc sampler
##'
##' This is a basic Metropolis-Hastings MCMC sampler.  The
##' \code{filter} is run with a set of parameters to evaluate the
##' likelihood. A new set of parameters is proposed, and these
##' likelihoods are compared, jumping with probability equal to their
##' ratio. This is repeated for \code{n_mcmc} proposals,
##' \code{n_chains} independent times.
##'
##' While this function is called \code{pmcmc} and requires a particle
##' filter object, there's nothing special about it for particle
##' filtering. However, we may need to add things in the future that
##' make assumptions about the particle filter, so we have named it
##' with a "p".
##'
##' @title Run a pmcmc sampler
##'
##' @param mcmc_range \code{data.frame} detailing parameters to sample.
##' Must contain columns 'name' (parameter names), 'init' (initial values),
##' 'min' (minimum values), 'max' (maximum values),
##' 'discrete (boolean indicating whether a discrete quantity)' and
##' 'target'
##'
##' @param lprior_funcs functions to calculate log prior for each parameter.
##' A named list for each parameter listed in \code{pars_to_sample}.
##' Each value must be a function which takes named parameter vector as
##'   input, returns a single numeric which is log of the prior probability.
##'
##' @param filter A particle filter object
##'
##' @param n_steps Number of MCMC steps
##'
##' @param proposal_kernel named matrix of proposal covariance for parameters
##'
##' @param return_proposals Logical indicating whether proposed parameter
##' jumps should be returned along with results
##'
##' @param n_chains Number of chains to run. If greater than 1, then
##'   we return a \code{mcstate_pmcmc_list} object which contains
##'   output for each chain, and an estimate of the Gelman's rhat.
##'
##' @param force_multichain Return a multichain object even when
##'   \code{n_chains} is 1
##'
##' @return Either a \code{mcstate_pmcmc} object or
##'   \code{mcstate_pmcmc_list} object, depending on the value of
##'   \code{n_chains}.
##'
##' @export
##' @import coda
pmcmc <- function(pars, filter, n_steps, return_proposals = FALSE,
                  n_chains = 1, force_multichain = FALSE) {
  assert_is(pars, "pmcmc_parameters")
  assert_is(filter, "particle_filter")
  assert_scalar_logical(return_proposals)

  target <- filter$run
  if (force_multichain || n_chains > 1) {
    mcmc_multichain(pars, target, n_steps, n_chains, return_proposals)
  } else {
    mcmc(pars, target, n_steps, return_proposals)
  }
}


##' Combine multiple chains into a single (master) chain
##'
##' @title Combine multiple chains into a single (master) chain
##'
##' @param x An mcstate_pmcmc_list containing chains to combine
##'
##' @param burn_in an integer denoting the number of samples to discard
##' from each chain
##'
##' @export
pmcmc_combine_chains <- function(x, burn_in) {
  assert_is(x, "mcstate_pmcmc_list")
  assert_scalar_positive_integer(burn_in)
  len <- nrow(x$chains[[1L]]$results)
  if (burn_in >= len) {
    stop("burn_in must be less than the total chain length")
  }
  drop <- seq_len(burn_in)
  ## TODO: should this include the chain id?
  ## TODO: Can we eliminate this function?
  ret <- do.call(rbind, lapply(x$chains, function(el) el$results[-drop, ]))
  rownames(ret) <- NULL
  ret
}


mcmc <- function(pars, target, n_steps, return_proposals) {
  curr_pars <- pars$initial()
  curr_lprior <- pars$prior(curr_pars)
  curr_llik <- target(pars$model(curr_pars))
  curr_lpost <- curr_lprior + curr_llik

  history <- pmcmc_history(n_steps, curr_pars, curr_lprior, curr_llik,
                           curr_lpost)

  if (return_proposals) {
    proposals <- pmcmc_history(n_steps, curr_pars, curr_lprior, curr_llik,
                               curr_lpost)
  }

  for (i in seq_len(n_steps)) {
    prop_pars <- pars$propose(curr_pars)
    prop_lprior <- pars$prior(prop_pars)
    prop_llik <- target(pars$model(prop_pars))
    prop_lpost <- prop_lprior + prop_llik

    exp(prop_lpost)
    exp(curr_lpost)

    if (runif(1) < exp(prop_lpost - curr_lpost)) {
      curr_pars <- prop_pars
      curr_lprior <- prop_lprior
      curr_llik <- prop_llik
      curr_lpost <- prop_lpost
    }

    history$add(curr_pars, curr_lprior, curr_llik, curr_lpost)
    if (return_proposals) {
      proposals$add(prop_pars, prop_lprior, prop_llik, prop_lpost)
    }
  }

  ## TODO: this is all up for grabs. Things that feel worthwhile:
  ## * do we return coda objects?
  ## * what parameter metadata do we save here (names, ranges, etc)
  ## * do we split the probabilities from the parameter estimates?
  results <- history$as_data_frame()
  out <- list(results = results,
              pars = pars$names(),
              acceptance_rate = acceptance_rate(results),
              ess = effective_size(results))

  if (return_proposals) {
    out$proposals <- proposals$as_data_frame()
  }

  class(out) <- "mcstate_pmcmc"
  out
}


mcmc_multichain <- function(pars, target, n_steps, n_chains,
                            return_proposals) {
  ## This could be parallelised with furrr::future_pmap. However, for
  ## the chains to be independent on threads, dust will need to
  ## support a long_jump call so that each can be advanced, and the
  ## random number generators will end up decoupled.  However, we've
  ## had problems with future etc failing to spawn properly, so
  ## holding off on that for now.
  chains <- lapply(seq_len(n_chains), function(i)
    mcmc(pars, target, n_steps, return_proposals))
  res <- list(rhat = gelman_diagnostic(chains, pars),
              chains = chains,
              pars = pars$summary())
  class(res) <- "mcstate_pmcmc_list"
  res
}


## This would be more nicely done as a simple R6 class but it's a bit
## slow in testing; this version speeds up the total mcmc runtime by a
## factor of ~3x (0.4s/1000 iterations to 0.13s/1000) mostly by
## reducing the number of garbage collections considerably.
pmcmc_history <- function(n, pars, log_prior, log_likelihood, log_posterior) {
  data <- vector("list", length = n + 1L)
  i <- 0L
  nms <- c(names(pars), "log_prior", "log_likelihood", "log_posterior")

  add <- function(pars, log_prior, log_likelihood, log_posterior) {
    i <<- i + 1L
    data[[i]] <<- c(pars, log_prior, log_likelihood, log_posterior)
  }

  as_data_frame <- function() {
    m <- matrix(unlist(data), ncol = length(nms), byrow = TRUE,
                dimnames = list(NULL, nms))
    as.data.frame(m, check.names = FALSE)
  }

  ## Add initial conditions
  add(pars, log_prior, log_likelihood, log_posterior)

  ## The rest is done in the pmcmc function
  list(add = add, as_data_frame = as_data_frame)
}


## Calculates the gelman diagnostic for multiple chains, if possible
gelman_diagnostic <- function(chains, pars) {
  chains_coda <- lapply(chains, function(x)
    coda::as.mcmc(x$results[pars$names()]))
  tryCatch(
    coda::gelman.diag(chains_coda),
    error = function(e) {
      message("Could not calculate rhat: ", e$message)
      NULL
    })
}


acceptance_rate <- function(chain) {
  ## TODO: this is actually pretty awful internally
  1 - coda::rejectionRate(coda::as.mcmc(chain))
}


effective_size <- function(chain) {
  ## TODO: do we ever want the ess of the probabilities?
  coda::effectiveSize(coda::as.mcmc(chain))
}
