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
##' @param n_chains Number of chains to run. If greater than 1, then
##'   we return a \code{mcstate_pmcmc_list} object which contains
##'   output for each chain, and an estimate of the Gelman's rhat.
##'
##' @param force_multichain Return a multichain object even when
##'   \code{n_chains} is 1
##'
##' @param collect_burnin Number of iterations to run before starting
##'   collection
##'
##' @param collect_interval Interval between samples
##'
##' @return Either a \code{mcstate_pmcmc} object or
##'   \code{mcstate_pmcmc_list} object, depending on the value of
##'   \code{n_chains}.
##'
##' @export
##' @import coda
pmcmc <- function(pars, filter, n_steps, return_state = FALSE, n_chains = 1,
                  force_multichain = FALSE) {
  assert_is(pars, "pmcmc_parameters")
  assert_is(filter, "particle_filter")
  assert_scalar_logical(return_proposals)

  reporting <- list(return_proposals = FALSE,
                    collect = list(
                      burnin = collect_burnin,
                      interval = collect_interval))

  target <- filter$run
  if (force_multichain || n_chains > 1) {
    mcmc_multichain(pars, target, n_steps, n_chains, reporting)
  } else {
    mcmc(pars, target, n_steps, reporting)
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


sample_trajectories <- function(history, index) {
  d <- dim(history)
  history <- history[, index, , drop = TRUE]
  dim(history) <- d[c(1, 3)]
  history
}


mcmc <- function(pars, filter, n_steps, save_state, save_trajectories) {
  history_pars <- history_collector(n_steps)
  history_probabilities <- history_collector(n_steps)
  history_state <- history_collector(n_steps)
  history_trajectories <- history_collector(n_steps)

  curr_pars <- pars$initial()
  curr_lprior <- pars$prior(curr_pars)
  curr_llik <- filter$run(pars$model(curr_pars), save_trajectories)
  curr_lpost <- curr_lprior + curr_llik

  history_pars$add(1L, curr_pars)
  history_probabilities$add(1L, c(curr_lprior, curr_llik, curr_lpost))

  n_particles <- filter$n_particles

  particle_idx <- sample.int(n_particles, 1)
  if (save_trajectories) {
    curr_trajectories <- sample_trajectories(filter$history(particle_idx))
    history_trajectories$add(1L, curr_trajectories)
  }

  if (save_state) {
    curr_state <- filter$state()[, particle_idx, drop = TRUE]
    history_state$add(1L, curr_state)
  }

  for (i in seq_len(n_steps)) {
    prop_pars <- pars$propose(curr_pars)
    prop_lprior <- pars$prior(prop_pars)
    prop_llik <- filter$run(pars$model(prop_pars), save_trajectories)
    prop_lpost <- prop_lprior + prop_llik

    if (runif(1) < exp(prop_lpost - curr_lpost)) {
      curr_pars <- prop_pars
      curr_lprior <- prop_lprior
      curr_llik <- prop_llik
      curr_lpost <- prop_lpost

      particle_idx <- sample.int(n_particles, 1)
      if (save_trajectories) {
        curr_trajectories <- sample_trajectories(filter$history(particle_idx))
      }
      if (save_state) {
        curr_state <- filter$state()[, particle_idx, drop = TRUE]
      }
    }

    idx <- i + 1L
    history_pars$add(idx, curr_pars)
    history_probabilities$add(idx, c(curr_lprior, curr_llik, curr_lpost))
    if (save_trajectories) {
      history_trajectories$add(idx, curr_trajectories)
    }
    if (save_state) {
      history_state$add(idx, curr_state)
    }
  }

  pars <- set_colnames(list_to_matrix(history_pars$get()),
                       names(curr_pars))
  probabilities <- set_colnames(list_to_matrix(history_probabilities$get()),
                       c("log_prior", "log_likelihood", "log_posterior"))

  state <- trajectories <- NULL
  if (save_state) {
    state <- t(list_to_matrix(history_state$get()))
  }
  if (save_trajectories) {
    trajectories <- list_to_array(history_trajectories$get())
  }

  out <- list(pars = pars, probabilities = probabilities, state = state,
              trajectories = trajectories)
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


## Generic history collector, collects anything at all into a list
history_collector <- function(n, burn_in = 0L, interval = 1L) {
  data <- vector("list", n + 1L)
  add <- function(i, value) {
    if (i > burn_in && i %% interval == 0) {
      data[[i]] <<- value
    }
  }

  get <- function() {
    data
  }

  list(add = add, get = get)
}


pmcmc_history <- function(n, reporting) {
  list(chain = history_collector(n_steps),
       state = history_collector(n_steps, reporting))
}
