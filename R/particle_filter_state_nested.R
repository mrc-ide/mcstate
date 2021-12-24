##' @title Nested particle filter state
##'
##' @description Nested particle filter internal state. This object is not
##'   ordinarily constructed directly by users, but via the
##'   `$run_begin` method to [mcstate::particle_filter]. It provides
##'   an advanced interface to the particle filter that allows
##'   partially running the particle filter over part of the time
##'   trajectory.
##'
##' This state object has a number of public fields that you can read
##'   but must not write (they are not read-only so you *could* write
##'   them, but don't).
particle_filter_state_nested <- R6::R6Class(
  "particle_filter_state_nested",
  cloneable = FALSE,

  private = list(
    generator = NULL,
    pars = NULL,
    data = NULL,
    data_split = NULL,
    steps = NULL,
    n_steps = NULL,
    n_particles = NULL,
    n_threads = NULL,
    initial = NULL,
    index = NULL,
    compare = NULL,
    gpu = NULL,
    save_restart_step = NULL,
    save_restart = NULL,
    min_log_likelihood = NULL
  ),

  public = list(
    ##' @field model The dust model generator being simulated (cannot be
    ##'   re-bound)
    model = NULL,

    ##' @field history The particle history, if created with
    ##'   `save_history = TRUE`. This is an internal format subject to
    ##    change.
    history = NULL,

    ##' @field restart_state Full model state at a series of points in
    ##'   time, if the model was created with non-`NULL` `save_restart`.
    ##'   This is a 4d array as described in [mcstate::particle_filter]
    restart_state = NULL,

    ##' @field log_likelihood The log-likelihood so far. This starts at
    ##'   0 when initialised and accumulates value for each step taken.
    log_likelihood = NULL,

    ##' @field log_likelihood_step The log-likelihood attributable to the
    ##' last step (i.e., the contribution to `log_likelihood` made on the
    ##' last call to `$step()`.
    log_likelihood_step = NULL,

    ##' @field current_step_index The index of the last completed step.
    current_step_index = 0L,

    ##' @description Initialise the particle filter state. Ordinarily
    ##' this should not be called by users, and so arguments are barely
    ##' documented.
    ##' @param pars Parameters for a single phase
    ##' @param generator A dust generator object
    ##' @param model If the generator has previously been initialised
    ##' @param data A [mcstate::particle_filter_data] data object
    ##' @param data_split The same data as `data` but split by step
    ##' @param steps A matrix of step beginning and ends
    ##' @param n_particles Number of particles to use
    ##' @param n_threads The number of threads to use
    ##' @param initial Initial condition function (or `NULL`)
    ##' @param index Index function (or `NULL`)
    ##' @param compare Compare function
    ##' @param gpu_config GPU configuration, passed to `generator`
    ##' @param seed Initial RNG seed
    ##' @param min_log_likelihood Early termination control
    ##' @param save_history Logical, indicating if we should save history
    ##' @param save_restart Vector of steps to save restart at
    initialize = function(pars, generator, model, data, data_split, steps,
                          n_particles, n_threads, initial, index, compare,
                          gpu_config, seed, min_log_likelihood,
                          save_history, save_restart) {
      ## These checks must exist differently for the two versions (and
      ## for the deterministic version)

      ## NOTE: this will generate a warning when updating docs but
      ## that's ok; see https://github.com/r-lib/roxygen2/issues/1067
      if (length(pars) < length(unique(data$population))) {
        stop(sprintf("'pars' must be at least the length of
                     'data$population', %d",
                     length(unique(data$population))))
      }

      if (is.null(model)) {
        pars_multi <- TRUE # DIFFERENCE
        model <- generator$new(pars = pars, step = steps[[1L]],
                               n_particles = n_particles, n_threads = n_threads,
                               seed = seed, gpu_config = gpu_config,
                               pars_multi = pars_multi)
        if (is.null(compare)) {
          model$set_index(integer(0))
          model$set_data(data_split)
        }
      } else {
        model$update_state(pars = pars, step = steps[[1L]])
      }

      if (!is.null(initial)) {
        initial_data <- pfsn_initial(model, initial, pars, n_particles)
        model$update_state(state = initial_data)
      }

      index_data <- pfsn_index(model, index)
      if (!is.null(compare) && !is.null(index_data$run)) {
        model$set_index(index_data$run)
      }

      if (save_history) {
        len <- nrow(steps) + 1L
        state <- model$state(index_data$state)
        history_value <- array(NA_real_, c(dim(state), len))
        array_last_dimension(history_value, 1) <- state
        history_order <- array(seq_len(n_particles),
                               c(model$shape(), len))
        self$history <- list(
          value = history_value,
          order = history_order,
          index = index_data$state)
      } else {
        self$history <- NULL
      }

      save_restart_step <- check_save_restart(save_restart, data)
      if (length(save_restart_step) > 0) {
        self$restart_state <- array(NA_real_,
                                    c(model$n_state(),
                                      model$shape(),
                                      length(save_restart)))
      } else {
        self$restart_state <- NULL
      }

      ## Constants
      private$generator <- generator
      private$pars <- pars
      private$data <- data
      private$data_split <- data_split
      private$steps <- steps
      private$n_particles <- n_particles
      private$n_threads <- n_threads
      private$initial <- initial
      private$index <- index
      private$compare <- compare
      private$gpu <- !is.null(gpu_config)
      private$min_log_likelihood <- min_log_likelihood
      private$save_restart_step <- save_restart_step
      private$save_restart <- save_restart

      ## Variable (see also history)
      self$model <- model
      self$log_likelihood <- rep(0.0, length(pars))
    },

    ##' @description Run the particle filter to the end of the data. This is
    ##' a convenience function around `$step()` which provides the correct
    ##' value of `step_index`
    run = function() {
      if (is.null(private$compare)) {
        particle_filter_compiled(self, private)
      } else {
        self$step(nrow(private$steps))
      }
    },

    ##' @description Take a step with the particle filter. This moves
    ##' the particle filter forward one step within the *data* (which
    ##' may correspond to more than one step with your model) and
    ##' returns the likelihood so far.
    ##'
    ##' @param step_index The step *index* to move to. This is not the same
    ##' as the model step, nor time, so be careful (it's the index within
    ##' the data provided to the filter). It is an error to provide
    ##' a value here that is lower than the current step index, or past
    ##' the end of the data.
    step = function(step_index) {
      steps <- private$steps
      curr <- self$current_step_index
      check_step(curr, step_index, private$steps, "Particle filter")

      model <- self$model
      compare <- private$compare

      ## This needs a little work in dust:
      ## https://github.com/mrc-ide/dust/issues/177
      if (is.null(compare)) {
        stop("Can't use low-level step with compiled particle filter (yet)")
      }

      steps <- private$steps
      data_split <- private$data_split
      pars <- private$pars

      restart_state <- self$restart_state
      save_restart_step <- private$save_restart_step
      save_restart <- !is.null(restart_state)

      history <- self$history
      save_history <- !is.null(history)
      save_history_index <- self$history$index
      history_value <- history$value
      history_order <- history$order

      log_likelihood <- self$log_likelihood
      n_particles <- private$n_particles

      min_log_likelihood <- private$min_log_likelihood

      for (t in seq(curr + 1L, step_index)) {
        step_end <- steps[t, 2L]
        state <- model$run(step_end)

        if (save_history) {
          array_last_dimension(history_value, t + 1L) <-
            model$state(save_history_index)
        }

        log_weights <- pfsn_compare(state, compare, data_split[[t]], pars)

        if (is.null(log_weights)) {
          if (save_history) {
            array_last_dimension(history_order, t + 1L) <- seq_len(n_particles)
          }
          log_likelihood_step <- NA_real_
        } else {
          tmp <- pfsn_weights(log_weights)
          log_likelihood_step <- tmp$log_likelihood

          log_likelihood <- log_likelihood + log_likelihood_step
          if (pfs_early_exit(log_likelihood, min_log_likelihood)) {
            log_likelihood <- -Inf
            break
          }

          kappa <- tmp$kappa
          model$reorder(kappa)
          if (save_history) {
            array_last_dimension(history_order, t + 1L) <- kappa
          }
        }

        if (save_restart) {
          i_restart <- match(step_end, save_restart_step)
          if (!is.na(i_restart)) {
            array_last_dimension(restart_state, i_restart) <- model$state()
          }
        }
      }

      ## Identical below here
      self$log_likelihood_step <- log_likelihood_step
      self$log_likelihood <- log_likelihood
      self$current_step_index <- step_index
      if (save_history) {
        history$value <- history_value
        history$order <- history_order
        self$history <- history
      }
      if (save_restart) {
        self$restart_state <- restart_state
      }

      log_likelihood
    },

    ##' @description Create a new `particle_filter_state_nested` object based
    ##' on this one (same model, position in time within the data) but
    ##' with new parameters. To do this, we create a new
    ##' `particle_filter_state_nested` with new parameters at the beginning of
    ##' the simulation (corresponding to the start of your data or the
    ##' `initial` argument to [mcstate::particle_filter]) with your new
    ##' `pars`, and then run the filter foward in time until it reaches
    ##' the same step as the parent model.
    ##'
    ##' @param pars New model parameters
    fork_smc2 = function(pars) {
      ## Identical except for the generator; that might go away?
      stopifnot(!private$gpu) # this won't work
      gpu_config <- NULL
      model <- NULL
      seed <- self$model$rng_state()
      save_history <- !is.null(self$history)

      ## Different generator used (pfsn not pfs)
      ret <- particle_filter_state_nested$new(
        pars, private$generator, model, private$data, private$data_split,
        private$steps, private$n_particles, private$n_threads,
        private$initial, private$index, private$compare, gpu_config,
        seed, private$min_log_likelihood, save_history, private$save_restart)

      ## Run it up to the same point
      ret$step(self$current_step_index)

      ## Set the seed in the parent model
      self$model$set_rng_state(ret$model$rng_state())

      ## Now, the user can use this model
      ret
    }
  ))


pfsn_initial <- function(model, initial, pars, n_particles) {
  state <- Map(initial, model$info(), rep(n_particles, length(pars)), pars)
  if (all(vlapply(state, is.null))) {
    return()
  }

  if (any(vlapply(state, is.list))) {
    stop("Setting 'step' from initial no longer supported")
  }

  len <- lengths(state)
  if (length(unique(len)) != 1) {
    stop(sprintf("initial() produced unequal state lengths %s",
                 str_collapse(len)))
  }

  if (is.null(dim(state[[1]]))) {
    state_array <- matrix(unlist(state, FALSE, FALSE), ncol = length(pars))
  } else {
    state_array <- list_to_array(state)
  }

  state_array
}


pfsn_index <- function(model, index) {
  if (is.null(index)) {
    return(NULL)
  }
  index_data <- lapply(model$info(), index)

  nok <- !all(vlapply(index_data[-1], identical, index_data[[1]]))
  if (nok) {
    stop("index must be identical across populations")
  }

  index_data[[1]]
}


## Can split and merge below here to come up with something better
## really.  All the wrangling here should come through so that compare
## does the entire thing through from comparison to weights,
## likelihood, but not kappa
pfsn_compare <- function(state, compare, data, pars) {
  log_weights <- lapply(seq_len(nlayer(state)), function(i)
    compare(array_drop(state[, , i, drop = FALSE], 3), data[[i]], pars[[i]]))

  if (all(lengths(log_weights) == 0)) {
    return(NULL)
  }
  vapply(log_weights, identity, numeric(length(log_weights[[1]])))
}


pfsn_weights <- function(log_weights) {
  n <- ncol(log_weights)
  weights <- lapply(seq_len(n), function(i)
    scale_log_weights(log_weights[, i]))
  log_likelihood <- vnapply(weights, "[[", "average")

  if (all(log_likelihood > -Inf)) {
    kappa <- vapply(weights, function(x) particle_resample(x$weights),
                    numeric(nrow(log_weights)))
  } else {
    ## No reordering
    kappa <- array(seq_len(nrow(log_weights)), dim(log_weights))
  }

  list(log_likelihood = log_likelihood, kappa = kappa)
}


pfsn_early_exit <- function(log_likelihood, min_log_likelihood) {
  scalar <- length(min_log_likelihood) == 1
  log_likelihood == -Inf ||
    ( scalar && sum(log_likelihood) < min_log_likelihood) ||
    (!scalar && all(log_likelihood < min_log_likelihood))
}
