##' Run a pmcmc sampler
##'
##' This is a basic Metropolis-Hastings MCMC sampler.  The
##' `filter` is run with a set of parameters to evaluate the
##' likelihood. A new set of parameters is proposed, and these
##' likelihoods are compared, jumping with probability equal to their
##' ratio. This is repeated for `n_mcmc` proposals.
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
##' @param n_steps Deprecated: use [mcstate::pmcmc_control] instead
##'
##' @param save_state Deprecated: use [mcstate::pmcmc_control] instead
##'
##' @param save_trajectories Deprecated: use [mcstate::pmcmc_control] instead
##'
##' @param progress Deprecated: use [mcstate::pmcmc_control] instead
##'
##' @param n_chains Deprecated: use [mcstate::pmcmc_control] instead
##'
##' @param initial Optional initial starting point. If given, it must
##'   be compatible with the parameters given in `pars`, and must be
##'   valid against your prior. You can use this to override the
##'   initial conditions saved in your `pars` object. You can provide
##'   either a vector of initial conditions, or a matrix with
##'   `n_chains` columns to use a different starting point for each
##'   chain.
##'
##' @param rerun_every Deprecated: use [mcstate::pmcmc_control] instead
##'
##' @param control A [mcstate::pmcmc_control] object which will set
##'   parameters. This will become the primary way of specifying
##'   options very soon and it is an error if you use both `control`
##'   and any of the parameters above aside from `pars` and `filter`.
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
pmcmc <- function(pars, filter, n_steps, save_state = TRUE,
                  save_trajectories = FALSE, progress = FALSE,
                  n_chains = 1, initial = NULL, rerun_every = Inf,
                  control = NULL) {

  assert_is(pars, c("pmcmc_parameters", "pmcmc_parameters_nested"))
  assert_is(filter, "particle_filter")

  if (is.null(control)) {
    warning("Please update your code to use pmcmc::pmcmc_control()",
            immediate. = TRUE)
    control <- pmcmc_control(n_steps,
                             n_chains = n_chains,
                             rerun_every = rerun_every,
                             save_state = save_state,
                             save_trajectories = save_trajectories,
                             progress = progress)
  } else {
    assert_is(control, "pmcmc_control")
    ok <- missing(n_steps) && missing(save_state) &&
      missing(save_trajectories) && missing(progress) && missing(n_chains) &&
      missing(rerun_every)
    if (!ok) {
      stop("Do not use deprecated arguments duplicated in pmcmc_control")
    }
  }

  if (inherits(pars, "pmcmc_parameters_nested")) {
    initial <- pmcmc_check_initial_nested(initial, pars, control$n_chains)
  } else {
    initial <- pmcmc_check_initial(initial, pars, control$n_chains)
  }


  if (control$n_workers == 1) {
    pmcmc_multiple_series(pars, initial, filter, control)
  } else {
    pmcmc_multiple_parallel(pars, initial, filter, control)
  }
}


pmcmc_single_chain <- function(pars, initial, filter, control, seed = NULL) {
  if (!is.null(seed)) {
    ## This will be triggered where control$use_parallel_seed is TRUE
    set.seed(seed$r)
    filter <- particle_filter_from_inputs(filter$inputs(), seed$dust)
  }
  obj <- pmcmc_state$new(pars, initial, filter, control)
  obj$run()
  obj$finish()
}

pmcmc_single_chain_nested <- function(pars, initial, filter, control,
                                      seed = NULL) {
  if (!is.null(seed)) {
    filter <- particle_filter_from_inputs(filter$inputs(), seed$dust)
  }
  obj <- pmcmc_state$new(pars, initial, filter, control)
  obj$run_nested()
  obj$finish_nested()
}


pmcmc_multiple_series <- function(pars, initial, filter, control) {
  if (control$use_parallel_seed) {
    seed <- make_seeds(control$n_chains, filter$inputs()$seed)
  } else {
    seed <- NULL
  }
  if (!is.null(control$n_threads_total)) {
    filter$set_n_threads(control$n_threads_total)
  }
  samples <- vector("list", control$n_chains)

  for (i in seq_along(samples)) {
    if (control$progress) {
      message(sprintf("Running chain %d / %d", i, control$n_chains))
    }
    if (inherits(pars, "pmcmc_parameters_nested")) {
      samples[[i]] <- pmcmc_single_chain_nested(pars, initial[, , i], filter,
                                                control, seed[[i]])
    } else {
      samples[[i]] <- pmcmc_single_chain(pars, initial[, i], filter, control,
                                         seed[[i]])
    }
  }
  if (length(samples) == 1) {
    samples[[1L]]
  } else {
    pmcmc_combine(samples = samples)
  }
}


pmcmc_multiple_parallel <- function(pars, initial, filter, control) {
  obj <- pmcmc_orchestrator$new(pars, initial, filter, control)
  obj$run()
  obj$finish()
}


## TODO: This does not check that the parameters are in range, or that
## they are appropriately discrete. We should add that in too at some
## point, though this overlaps with some outstanding validation in the
## smc2 branch.
pmcmc_check_initial <- function(initial, pars, n_chains) {
  nms <- pars$names()
  n_pars <- length(nms)
  if (is.null(initial)) {
    initial <- pars$initial()
  }
  if (is.matrix(initial)) {
    if (nrow(initial) != n_pars) {
      stop(sprintf("Expected a matrix with %d rows for 'initial'", n_pars))
    }
    if (ncol(initial) != n_chains) {
      stop(sprintf("Expected a matrix with %d columns for 'initial'", n_chains))
    }
    if (!is.null(rownames(initial)) && !identical(rownames(initial), nms)) {
      stop("If 'initial' has rownames, they must match pars$names()")
    }
    ok <- apply(initial, 2, function(p) is.finite(pars$prior(p)))
    if (any(!ok)) {
      stop(sprintf(
        "Starting point does not have finite prior probability (%s)",
        paste(which(!ok), collapse = ", ")))
    }
  } else {
    if (length(initial) != n_pars) {
      stop(sprintf("Expected a vector of length %d for 'initial'", n_pars))
    }
    if (!is.null(names(initial)) && !identical(names(initial), nms)) {
      stop("If 'initial' has names, they must match pars$names()")
    }
    if (!is.finite(pars$prior(initial))) {
      stop("Starting point does not have finite prior probability")
    }
    initial <- matrix(initial, n_pars, n_chains)
  }
  dimnames(initial) <- list(nms, NULL)
  initial
}

## TODO: This does not check that the parameters are in range, or that
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
    if (nlayer(initial) != n_chains) {
      stop(sprintf("Expected an array with %d layers for 'initial'", n_chains))
    }
    if (ncol(initial) != n_pars) {
      stop(sprintf("Expected an array with %d columns for 'initial'", n_pars))
    }
    if (nrow(initial) != n_pops) {
      stop(sprintf("Expected an array with %d rows for 'initial'", n_pops))
    }
    if (!is.null(rownames(initial)) && !identical(rownames(initial), pops)) {
      stop("If 'initial' has rownames, they must match pars$populations()")
    }
    if (!is.null(colnames(initial)) && !identical(colnames(initial), nms)) {
      stop("If 'initial' has colnames, they must match pars$names()")
    }

    dimnames(initial) <- list(pops, nms, NULL)

    ok <- apply(initial, 3, function(p) all(is.finite(pars$prior(p))))
    if (any(!ok)) {
      stop(sprintf(
        "Starting point does not have finite prior probability (%s)",
        paste(which(!ok), collapse = ", ")))
    }
  } else {
    if (NCOL(initial) != n_pars) {
      stop(sprintf("Expected a matrix with %d columns for 'initial'", n_pars))
    }
    if (NROW(initial) != n_pops) {
      stop(sprintf("Expected a matrix with %d rows for 'initial'", n_pops))
    }
    if (!is.null(rownames(initial)) && !identical(rownames(initial), pops)) {
      stop("If 'initial' has rownames, they must match pars$populations()")
    }
    if (!is.null(colnames(initial)) && !identical(colnames(initial), nms)) {
      stop("If 'initial' has colnames, they must match pars$names()")
    }

    dimnames(initial) <- list(pops, nms)

    if (any(!is.finite(pars$prior(initial)))) {
      stop("Starting point does not have finite prior probability")
    }

    initial <- array(initial, c(n_pops, n_pars, n_chains),
                     dimnames = list(pops, nms, NULL))
  }

  initial
}
