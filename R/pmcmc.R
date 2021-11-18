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
##' @param initial Optional initial starting point. If given, it must
##'   be compatible with the parameters given in `pars`, and must be
##'   valid against your prior. You can use this to override the
##'   initial conditions saved in your `pars` object. You can provide
##'   either a vector of initial conditions, or a matrix with
##'   `n_chains` columns to use a different starting point for each
##'   chain.
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
pmcmc <- function(pars, filter, initial = NULL, control = NULL) {
  assert_is(pars, c("pmcmc_parameters", "pmcmc_parameters_nested"))
  assert_is(filter, c("particle_filter", "particle_deterministic"))
  assert_is(control, "pmcmc_control")

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


##' Run a pMCMC, with sensible random number behaviour, but schedule
##' execution of the chains yourself. Use this if you want to
##' distribute chains over (say) the nodes of an HPC system.
##'
##' Basic usage will look like
##'
##' ```
##' inputs <- mcstate::pmcmc_chains_prepare(pars, filter, control = control)
##' samples_data <-
##'   lapply(seq_len(control$n_chains), mcstate::pmcmc_chains_run, inputs)
##' samples <- mcstate::pmcmc_combine(samples = samples)
##' ```
##'
##' You can safely parallelise (or not) however you like at the
##' `lapply` call and get the same outputs regardless.
##'
##' @title pMCMC with manual chain scheduling
##'
##' @inheritParams pmcmc
##' @export
pmcmc_chains_prepare <- function(pars, filter, initial = NULL, control = NULL) {
  assert_is(pars, c("pmcmc_parameters", "pmcmc_parameters_nested"))
  assert_is(filter, "particle_filter")
  assert_is(control, "pmcmc_control")

  if (control$n_workers != 1) {
    stop("'n_workers' must be 1")
  }
  if (!control$use_parallel_seed) {
    stop("'use_parallel_seed' must be TRUE")
  }

  if (inherits(pars, "pmcmc_parameters_nested")) {
    initial <- pmcmc_check_initial_nested(initial, pars, control$n_chains)
  } else {
    initial <- pmcmc_check_initial(initial, pars, control$n_chains)
  }

  seed <- make_seeds(control$n_chains, filter$inputs()$seed, filter$model)

  ret <- list(pars = pars, initial = initial, filter = filter,
              control = control, seed = seed)
  class(ret) <- "pmcmc_inputs"
  ret
}


##' @param chain_id The chain index to run (1, 2, ..., `control$n_chains`)
##'
##' @param inputs A `pmcmc_inputs` object created by `pmcmc_chains_prepare`
##'
##' @param path Optionally a directory to save output in. This might
##'   be useful if splitting work across multiple processes. Samples
##'   will be saved at `<path>/samples_<chain_id>.rds`
##'
##' @export
##' @rdname pmcmc_chains_prepare
pmcmc_chains_run <- function(chain_id, inputs, path = NULL) {
  assert_is(inputs, "pmcmc_inputs")
  assert_scalar_positive_integer(chain_id)
  if (chain_id < 1 || chain_id > inputs$control$n_chains) {
    stop(sprintf("'chain_id' must be an integer in 1..%d",
                 inputs$control$n_chains))
  }
  samples <- pmcmc_run_chain(chain_id, inputs$pars, inputs$initial,
                             inputs$filter, inputs$control, inputs$seed)
  if (is.null(path)) {
    samples
  } else {
    fmt <- "samples_%d.rds"
    dir.create(path, FALSE, TRUE)
    dest <- file.path(path, sprintf(fmt, chain_id))
    saveRDS(samples, dest)
    dest
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
    set.seed(seed$r)
    filter <- particle_filter_from_inputs(filter$inputs(), seed$dust)
  }
  obj <- pmcmc_state$new(pars, initial, filter, control)
  obj$run_nested()
  obj$finish_nested()
}


pmcmc_multiple_series <- function(pars, initial, filter, control) {
  if (control$use_parallel_seed) {
    seed <- make_seeds(control$n_chains, filter$inputs()$seed, filter$model)
  } else {
    seed <- NULL
  }
  samples <- vector("list", control$n_chains)

  ## Ensure that even if the control object has been updated we don't
  ## do anything impossible:
  control$n_steps_each <- control$n_steps

  for (i in seq_along(samples)) {
    samples[[i]] <- pmcmc_run_chain(i, pars, initial, filter, control, seed)
  }

  if (length(samples) == 1) {
    samples[[1L]]
  } else {
    pmcmc_combine(samples = samples)
  }
}


pmcmc_run_chain <- function(chain_id, pars, initial, filter, control, seed) {
  if (control$progress) {
    message(sprintf("Running chain %d / %d", chain_id, control$n_chains))
  }
  if (!is.null(control$n_threads_total)) {
    filter$set_n_threads(control$n_threads_total)
  }
  if (inherits(pars, "pmcmc_parameters_nested")) {
    pmcmc_single_chain_nested(pars, initial[, , chain_id], filter,
                              control, seed[[chain_id]])
  } else {
    pmcmc_single_chain(pars, initial[, chain_id], filter, control,
                       seed[[chain_id]])
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
