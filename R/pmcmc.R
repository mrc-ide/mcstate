##' Run a pmcmc sampler
##'
##' This is a basic Metropolis-Hastings MCMC sampler.  The
##' `filter` is run with a set of parameters to evaluate the
##' likelihood. A new set of parameters is proposed, and these
##' likelihoods are compared, jumping with probability equal to their
##' ratio. This is repeated for `n_steps` proposals.
##'
##' While this function is called `pmcmc` and requires a particle
##' filter object, there's nothing special about it for particle
##' filtering. However, we may need to add things in the future that
##' make assumptions about the particle filter, so we have named it
##' with a "p".
##'
##' @title Run a pmcmc sampler
##'
##' @param pars A [`pmcmc_parameters`] object containing
##'   information about parameters (ranges, priors, proposal kernel,
##'   translation functions for use with the particle filter).
##'
##' @param filter A [`particle_filter`] object
##'
##' @param initial Optional initial starting point. If given, it must
##'   be compatible with the parameters given in `pars`, and must be
##'   valid against your prior. You can use this to override the
##'   initial conditions saved in your `pars` object. You can provide
##'   either a vector of initial conditions, or a matrix with
##'   `n_chains` columns to use a different starting point for each
##'   chain.
##'
##' @param control A [mcstate::pmcmc_control] object which will
##'   control how the MCMC runs, including the number of steps etc.
##'
##' @return A `mcstate_pmcmc` object containing `pars`
##'   (sampled parameters) and `probabilities` (log prior, log
##'   likelihood and log posterior values for these
##'   probabilities). Two additional fields may be present:
##'   `state` (if `return_state` was `TRUE`),
##'   containing the final state of a randomly selected particle at
##'   the end of the simulation, for each step (will be a matrix with
##'   as many rows as your state has variables, and as `n_steps +
##'   1` columns corresponding to each step). `trajectories` will
##'   include a 3d array of particle trajectories through the
##'   simulation (if `return_trajectories` was `TRUE`).
##'
##' @export
pmcmc <- function(pars, filter, initial = NULL, control = NULL) {
  assert_is(pars, c("pmcmc_parameters", "pmcmc_parameters_nested"))
  assert_is(filter, c("particle_filter", "particle_deterministic"))
  assert_is(control, "pmcmc_control")
  pmcmc_check_control(control)

  if (control$n_workers == 1) {
    pmcmc_multiple_series(pars, initial, filter, control)
  } else {
    pmcmc_multiple_parallel(pars, initial, filter, control)
  }
}


pmcmc_multiple_series <- function(pars, initial, filter, control) {
  initial <- pmcmc_check_initial(initial, pars, control$n_chains)
  if (control$use_parallel_seed) {
    seed <- make_seeds(control$n_chains, filter$inputs()$seed, filter$model)
  } else {
    seed <- NULL
  }

  samples <- vector("list", control$n_chains)

  for (i in seq_along(samples)) {
    if (control$progress) {
      message(sprintf("Running chain %d / %d", i, control$n_chains))
    }
    if (control$use_parallel_seed) {
      ## TODO (#174): this feels pretty weird; can we just add a
      ## "set_rng_state" method for the filter? (see similar code in
      ## pmcmc_chains)
      set.seed(seed[[i]]$r)
      filter <- particle_filter_from_inputs(filter$inputs(), seed[[i]]$dust)
    }
    samples[[i]] <- pmcmc_run_chain(pars, initial[[i]], filter, control, NULL)
  }

  if (control$n_chains == 1) {
    samples[[1L]]
  } else {
    pmcmc_combine(samples = samples)
  }
}


pmcmc_run_chain <- function(pars, initial, filter, control, n_threads) {
  if (!is.null(n_threads)) {
    filter$set_n_threads(n_threads)
  } else if (!is.null(control$n_threads_total)) {
    filter$set_n_threads(control$n_threads_total / control$n_workers)
  }
  obj <- pmcmc_state$new(pars, initial, filter, control)
  obj$run()
  obj$finish()
}


pmcmc_multiple_parallel <- function(pars, initial, filter, control) {
  obj <- pmcmc_orchestrator$new(pars, initial, filter, control)
  obj$run()
  obj$finish()
}


pmcmc_check_initial <- function(initial, pars, n_chains) {
  if (inherits(pars, "pmcmc_parameters_nested")) {
    initial <- pmcmc_check_initial_nested(initial, pars, n_chains)
  } else {
    initial <- pmcmc_check_initial_simple(initial, pars, n_chains)
  }

  ## Process that into a list so that access later is simple (and
  ## identical regardless of nestedness or not)
  lapply(seq_len(n_chains), function(i)
    array_drop(array_last_dimension(initial, i), length(dim(initial))))
}


## TODO (#175): This does not check that the parameters are in range,
## or that they are appropriately discrete. We should add that in too
## at some point, though this overlaps with some outstanding
## validation in the smc2 branch.
pmcmc_check_initial_simple <- function(initial, pars, n_chains) {
  nms <- pars$names()
  n_pars <- length(nms)
  if (is.null(initial)) {
    initial <- pars$initial()
  }
  if (is.matrix(initial)) {
    assert_dimensions(initial, c(n_pars, n_chains))
    initial <- assert_dimnames(initial, list(parameters = nms, NULL))
  } else {
    assert_dimensions(initial, n_pars)
    assert_dimnames(initial, list(parameters = nms))
    initial <- array(initial, c(n_pars, n_chains),
                     dimnames = list(nms, NULL))
  }

  summary <- pars$summary()
  ok <- apply(initial, 2, function(p) all(p >= summary$min))
  if (any(!ok)) {
    stop(sprintf(
      "'initial' is less than 'min' (%s) (chain %s)",
      paste(summary$min, collapse = ", "),
      paste(which(!ok), collapse = ", ")))
  }
  ok <- apply(initial, 2, function(p) all(p <= summary$max))
  if (any(!ok)) {
    stop(sprintf(
      "'initial' is greater than 'max' (%s) (chain %s)",
      paste(summary$max, collapse = ", "),
      paste(which(!ok), collapse = ", ")))
  }
  ok <- apply(initial, 2, function(p) is.finite(pars$prior(p)))
  if (any(!ok)) {
    stop(sprintf(
      "Starting point does not have finite prior probability (chain %s)",
      paste(which(!ok), collapse = ", ")))
  }

  initial
}

## TODO (#175): This does not check that the parameters are in range, or that
## they are appropriately discrete. We should add that in too at some
## point, though this overlaps with some outstanding validation in the
## smc2 branch.
pmcmc_check_initial_nested <- function(initial, pars, n_chains) {
  nms <- pars$names()
  pops <- pars$populations()
  n_pars <- length(nms)
  n_pops <- length(pops)

  if (is.null(initial)) {
    initial <- pars$initial()
  }
  if (is_3d_array(initial)) {
    assert_dimensions(initial, c(n_pars, n_pops, n_chains))
    initial <- assert_dimnames(
      initial, list(parameters = nms, populations = pops, NULL))
  } else {
    assert_is(initial, "matrix")
    assert_dimensions(initial, c(n_pars, n_pops))
    assert_dimnames(initial, list(parameters = nms, populations = pops))
    initial <- array(initial, c(n_pars, n_pops, n_chains),
                     dimnames = list(nms, pops, NULL))
  }

  ok <- apply(initial, 3, function(p) all(is.finite(pars$prior(p))))
  if (any(!ok)) {
    stop(sprintf(
      "Starting point does not have finite prior probability (chain %s)",
      paste(which(!ok), collapse = ", ")))
  }

  initial
}
