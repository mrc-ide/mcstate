##' Run a pmcmc sampler
##'
##' This is a basic Metropolis-Hastings MCMC sampler.  The
##' \code{filter} is run with a set of parameters to evaluate the
##' likelihood. A new set of parameters is proposed, and these
##' likelihoods are compared, jumping with probability equal to their
##' ratio. This is repeated for \code{n_mcmc} proposals.
##'
##' While this function is called \code{pmcmc} and requires a particle
##' filter object, there's nothing special about it for particle
##' filtering. However, we may need to add things in the future that
##' make assumptions about the particle filter, so we have named it
##' with a "p".
##'
##' @title Run a pmcmc sampler
##'
##' @param pars A \code{\link{pmcmc_parameters}} object containing
##'   information about parameters (ranges, priors, proposal kernel,
##'   translation functions for use with the particle filter).
##'
##' @param filter A \code{\link{particle_filter}} object
##'
##' @param n_steps Number of MCMC steps to run
##'
##' @param save_state Logical, indicating if the state should be saved
##'   at the end of the simulation. If \code{TRUE}, then a single
##'   randomly selected particle's state will be collected at the end
##'   of each mcmc step. This is the full state (i.e., unaffected by
##'   and \code{index} used in the particle filter) so that the
##'   process may be restarted from this point for projections.  If
##'   \code{save_trajectories} is \code{TRUE} the same particle will
##'   be selected for each.
##'
##' @param save_trajectories Logical, indicating if the particle
##'   trajectories should be saved during the simulation. If
##'   \code{TRUE}, then a single randomly selected particle's
##'   trajectory will be collected at the end of each mcmc step.  This
##'   is the filtered state (i.e., using the \code{state} component of
##'   \code{index} provided to the particle filter).  If
##'   \code{save_state} is \code{TRUE} the same particle will
##'   be selected for each.
##'
##' @return A \code{mcstate_pmcmc} object containing \code{pars}
##'   (sampled parameters} and \code{probabilities} (log prior, log
##'   likelihood and log posterior values for these
##'   probabilities). Two additional fields may be present:
##'   \code{state} (if \code{return_state} was \code{TRUE}),
##'   containing the final state of a randomly selected particle at
##'   the end of the simulation, for each step (will be a matrix with
##'   as many rows as your state has variables, and as \code{n_steps +
##'   1} columns corresponding to each step). \code{trajectories} will
##'   include a 3d array of particle trajectories through the
##'   simulation (if \code{return_trajectories} was \code{TRUE}).
##'
##' @export
pmcmc <- function(pars, filter, n_steps, save_state = FALSE,
                  save_trajectories = FALSE) {
  assert_is(pars, "pmcmc_parameters")
  assert_is(filter, "particle_filter")
  assert_scalar_positive_integer(n_steps)
  assert_scalar_logical(save_state)
  assert_scalar_logical(save_trajectories)

  n_particles <- filter$n_particles

  history_pars <- history_collector(n_steps)
  history_probabilities <- history_collector(n_steps)
  history_state <- history_collector(n_steps)
  history_trajectories <- history_collector(n_steps)

  curr_pars <- pars$initial()
  curr_lprior <- pars$prior(curr_pars)
  curr_llik <- filter$run(pars$model(curr_pars), save_trajectories)
  curr_lpost <- curr_lprior + curr_llik

  history_pars$add(curr_pars)
  history_probabilities$add(c(curr_lprior, curr_llik, curr_lpost))

  particle_idx <- sample.int(n_particles, 1)

  if (save_trajectories) {
    curr_trajectories <- sample_trajectory(filter$history(particle_idx))
    history_trajectories$add(curr_trajectories)
  }

  if (save_state) {
    curr_state <- filter$state()[, particle_idx, drop = TRUE]
    history_state$add(curr_state)
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
        curr_trajectories <- sample_trajectory(filter$history(particle_idx))
      }
      if (save_state) {
        curr_state <- filter$state()[, particle_idx, drop = TRUE]
      }
    }

    history_pars$add(curr_pars)
    history_probabilities$add(c(curr_lprior, curr_llik, curr_lpost))
    if (save_trajectories) {
      history_trajectories$add(curr_trajectories)
    }
    if (save_state) {
      history_state$add(curr_state)
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

  mcstate_pmcmc(pars, probabilities, state, trajectories)
}


## A utility function for sampling a trajectory and safely dropping
## the dimensionality even if there is only one state vector
sample_trajectory <- function(history, index) {
  ret <- history[, index, , drop = TRUE]
  if (is.null(dim(ret))) {
    dim(ret) <- dim(history)[c(1, 3)]
  }
  ret
}


## Generic history collector, collects anything at all into a list
##
## This would be more nicely done as a simple R6 class but it's a bit
## slow in testing; this version speeds up the total mcmc runtime by a
## factor of ~3x (0.4s/1000 iterations to 0.13s/1000) mostly by
## reducing the number of garbage collections considerably.
history_collector <- function(n) {
  data <- vector("list", n + 1L)
  i <- 0L
  add <- function(value) {
    i <<- i + 1L
    data[[i]] <<- value
  }

  get <- function() {
    data
  }

  list(add = add, get = get)
}
