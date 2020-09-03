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
##' @param n_steps Number of MCMC steps to run
##'
##' @param save_state Logical, indicating if the state should be saved
##'   at the end of the simulation. If `TRUE`, then a single
##'   randomly selected particle's state will be collected at the end
##'   of each MCMC step. This is the full state (i.e., unaffected by
##'   and `index` used in the particle filter) so that the
##'   process may be restarted from this point for projections.  If
##'   `save_trajectories` is `TRUE` the same particle will
##'   be selected for each. The default is `TRUE`, which will
##'   cause `n_state` * `n_steps` of data to be output
##'   alongside your results. Set this argument to `FALSE` to
##'   save space, or use [pmcmc_thin()] after running the
##'   MCMC.
##'
##' @param save_trajectories Logical, indicating if the particle
##'   trajectories should be saved during the simulation. If
##'   `TRUE`, then a single randomly selected particle's
##'   trajectory will be collected at the end of each MCMC step.  This
##'   is the filtered state (i.e., using the `state` component of
##'   `index` provided to the particle filter).  If
##'   `save_state` is `TRUE` the same particle will
##'   be selected for each.
##'
##' @param progress Logical, indicating if a progress bar should be
##'   displayed, using [`progress::progress_bar`].
##'
##' @param n_chains Optional integer, indicating the number of chains
##'   to run. If more than one then we run a series of chains and
##'   merge them with [pmcmc_combine()]. Chains are run in series,
##'   with the same filter.
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
                  n_chains = 1) {
  assert_is(pars, "pmcmc_parameters")
  assert_is(filter, "particle_filter")
  assert_scalar_positive_integer(n_steps)
  assert_scalar_logical(save_state)
  assert_scalar_logical(save_trajectories)
  assert_scalar_positive_integer(n_chains)

  if (n_chains == 1) {
    pmcmc_single_chain(pars, filter, n_steps,
                       save_state, save_trajectories, progress)
  } else {
    samples <- vector("list", n_chains)
    for (i in seq_along(samples)) {
      if (progress) {
        message(sprintf("Running chain %d / %d", i, n_chains))
      }
      samples[[i]] <- pmcmc_single_chain(pars, filter, n_steps,
                                         save_state, save_trajectories,
                                         progress)
    }
    pmcmc_combine(samples = samples)
  }
}

pmcmc_single_chain <- function(pars, filter, n_steps,
                               save_state, save_trajectories,
                               progress) {
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

  tick <- pmcmc_progress(n_steps, progress)

  for (i in seq_len(n_steps)) {
    tick()

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

  pars_matrix <- set_colnames(list_to_matrix(history_pars$get()),
                              names(curr_pars))
  probabilities <- set_colnames(list_to_matrix(history_probabilities$get()),
                       c("log_prior", "log_likelihood", "log_posterior"))

  predict <- state <- trajectories <- NULL

  if (save_state) {
    state <- t(list_to_matrix(history_state$get()))

    ## Extract the transformation function from the pars object as
    ## it's all we want; this is a bit of a hack but it works for now
    ## at least.
    transform <- pars[[".__enclos_env__"]]$private$transform

    ## Extract compare function, data and particle filter index for use in
    ## re-running particle filter
    compare <- filter[[".__enclos_env__"]]$private$compare
    data <- filter[[".__enclos_env__"]]$private$data
    filter_index <- filter[[".__enclos_env__"]]$private$index

    ## Information about how the particle filter was configured:
    info <- filter$predict_info()

    ## This information is required in order to run predictions but
    ## adds a significant space overhead to the saved data (model is
    ## fairly large, and transform could capture anything in scope if
    ## it is a closure).
    predict <- list(transform = transform,
                    model = filter$model,
                    n_threads = info$n_threads,
                    index = info$index,
                    rate = info$rate,
                    seed = info$seed,
                    step = last(info$step),
                    compare = compare,
                    data = data, 
                    filter_index = filter_index)
  }

  if (save_trajectories) {
    info <- filter$predict_info()
    ## Permute trajectories from [state x mcmc x particle] to
    ## [state x particle x mcmc] so that they match the ones that we
    ## will generate with predict
    trajectories_state <-
      aperm(list_to_array(history_trajectories$get()), c(1, 3, 2))
    rownames(trajectories_state) <- names(info$index)
    trajectories <- mcstate_trajectories(info$step, info$rate,
                                         trajectories_state, FALSE)
  }

  mcstate_pmcmc(pars_matrix, probabilities, state, trajectories, predict)
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
