##' @title Particle filter state
##'
##' @description Particle filter internal state. This object is not
##'   ordinarily constructed directly by users, but via the
##'   `$run_begin` method to [mcstate::particle_filter]. It provides
##'   an advanced interface to the particle filter that allows
##'   partially running the particle filter over part of the time
##'   trajectory.
##'
##' This state object has a number of public fields that you can read
##'   but must not write (they are not read-only so you *could* write
##'   them, but don't).
particle_filter_state <- R6::R6Class(
  "particle_filter_state",
  cloneable = FALSE,

  private = list(
    generator = NULL,
    pars = NULL,
    data = NULL,
    data_split = NULL,
    steps = NULL,
    n_particles = NULL,
    n_threads = NULL,
    initial = NULL,
    index = NULL,
    compare = NULL,
    gpu = NULL,
    save_restart_step = NULL,
    save_restart = NULL,
    min_log_likelihood = NULL,
    support = NULL,

    run_compiled = function() {
      history <- self$history
      save_history <- !is.null(history)
      save_history_index <- self$history$index

      model <- self$model
      if (save_history) {
        model$set_index(save_history_index)
        on.exit(model$set_index(integer(0)))
      }

      res <- model$filter(save_history, private$save_restart_step)

      self$log_likelihood_step <- NA_real_
      self$log_likelihood <- res$log_likelihood
      self$current_step_index <- nrow(private$steps)
      if (save_history) {
        self$history <- list(value = res$trajectories,
                             index = save_history_index)
      }
      self$restart_state <- res$snapshots
      res$log_likelihood
    }
  ),

  public = list(
    ##' @field model The dust model being simulated
    model = NULL,

    ##' @field history The particle history, if created with
    ##'   `save_history = TRUE`. This is an internal format subject to
    ##    change.
    history = NULL,

    ##' @field restart_state Full model state at a series of points in
    ##'   time, if the model was created with non-`NULL` `save_restart`.
    ##'   This is a 3d (or greater) array as described in
    ##'   [mcstate::particle_filter]
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
      pars_multi <- inherits(data, "particle_filter_data_nested")
      support <- particle_filter_state_support(pars_multi)

      if (pars_multi) {
        if (!is.null(names(pars))) {
          stop("Expected an unnamed list of parameters")
        }
        if (length(pars) != attr(data, "n_populations")) {
          stop(sprintf("'pars' must have length %d (following data$%s)",
                       attr(data, "n_populations"), attr(data, "population")))
        }
      }

      if (is.null(model)) {
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
        initial_data <- support$initial(model, initial, pars, n_particles)
        model$update_state(state = initial_data)
      }

      index_data <- support$index(model, index)
      if (!is.null(compare) && !is.null(index_data$run)) {
        model$set_index(index_data$run)
      }

      ## The model shape is [n_particles, <any multi-par structure>]
      shape <- model$shape()

      if (save_history) {
        len <- nrow(steps) + 1L
        state <- model$state(index_data$state)
        history_value <- array(NA_real_, c(dim(state), len))
        array_last_dimension(history_value, 1) <- state
        history_order <- array(seq_len(n_particles), c(shape, len))
        self$history <- list(
          value = history_value,
          order = history_order,
          index = index_data$state)
      } else {
        self$history <- NULL
      }

      save_restart_step <- check_save_restart(save_restart, data)
      if (length(save_restart_step) > 0) {
        self$restart_state <-
          array(NA_real_, c(model$n_state(), shape, length(save_restart)))
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
      private$support <- support

      ## Variable (see also history)
      self$model <- model
      self$log_likelihood <- rep(0, prod(shape[-1]))
    },

    ##' @description Run the particle filter to the end of the data. This is
    ##' a convenience function around `$step()` which provides the correct
    ##' value of `step_index`
    run = function() {
      if (is.null(private$compare)) {
        private$run_compiled()
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
    ##'
    ##' @param partial Logical, indicating if we should return the partial
    ##' likelihood, due to this step, rather than the full likelihood so far.
    step = function(step_index, partial = FALSE) {
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
      support <- private$support

      for (t in seq(curr + 1L, step_index)) {
        step_end <- steps[t, 2L]
        state <- model$run(step_end)

        if (save_history) {
          array_last_dimension(history_value, t + 1L) <-
            model$state(save_history_index)
        }

        weights <- support$compare(state, compare, data_split[[t]], pars)

        if (is.null(weights)) {
          if (save_history) {
            array_last_dimension(history_order, t + 1L) <- seq_len(n_particles)
          }
          log_likelihood_step <- NA_real_
        } else {
          log_likelihood_step <- weights$average
          log_likelihood <- log_likelihood + log_likelihood_step

          if (particle_filter_early_exit(log_likelihood, min_log_likelihood)) {
            log_likelihood <- -Inf
            break
          }

          kappa <- particle_resample(weights$weights)
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

      if (partial) {
        log_likelihood_step
      } else {
        log_likelihood
      }
    },

    ##' @description Create a new `particle_filter_state` object based on
    ##' this one (same model, position in time within the data) but with
    ##' new parameters, to support the "multistage particle filter".
    ##' Unlike `fork_smc2`, here the parameters may imply a different
    ##' model shape and arbitrary transformations of the state are
    ##' allowed.  The model is not rerun to the current point, just
    ##' transformed at that point.
    ##'
    ##' @param pars New model parameters
    ##'
    ##' @param transform_state A function to transform the model state
    ##'   from the old to the new parameter set.  See
    ##'   [mcstate::multistage_epoch()] for details.
    fork_multistage = function(pars, transform_state) {
      stopifnot(!private$gpu) # this won't work
      gpu_config <- NULL
      model <- NULL
      seed <- self$model$rng_state()
      save_history <- !is.null(self$history)
      initial <- NULL

      if (is.null(pars)) {
        pars <- self$model$pars()
      }
      ret <- particle_filter_state$new(
        pars, private$generator, model, private$data, private$data_split,
        private$steps, private$n_particles, private$n_threads,
        initial, private$index, private$compare, gpu_config,
        seed, private$min_log_likelihood, save_history, private$save_restart)

      state <- transform_state(self$model$state(), self$model, ret$model)
      step <- self$model$step()

      ret$model$update_state(state = state, step = step)
      ret$current_step_index <- self$current_step_index
      ret$log_likelihood <- self$log_likelihood
      ret$log_likelihood_step <- self$log_likelihood_step

      ret
    },

    ##' @description Create a new `particle_filter_state` object based
    ##' on this one (same model, position in time within the data) but
    ##' with new parameters, run up to the date, to support the [smc2()]
    ##' algorithm. To do this, we create a new
    ##' `particle_filter_state` with new parameters at the beginning of
    ##' the simulation (corresponding to the start of your data or the
    ##' `initial` argument to [mcstate::particle_filter]) with your new
    ##' `pars`, and then run the filter foward in time until it reaches
    ##' the same step as the parent model.
    ##'
    ##' @param pars New model parameters
    fork_smc2 = function(pars) {
      stopifnot(!private$gpu) # this won't work
      gpu_config <- NULL
      model <- NULL
      seed <- self$model$rng_state()
      save_history <- !is.null(self$history)

      ret <- particle_filter_state$new(
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


## Used for both the normal and deterministic particle filter
check_step <- function(curr, step_index, steps, name) {
  n_steps <- nrow(steps)
  if (curr >= n_steps) {
    stop(sprintf("%s has reached the end of the data", name))
  }
  if (step_index > n_steps) {
    stop(sprintf("step_index %d is beyond the length of the data (max %d)",
                 step_index, n_steps))
  }
  if (step_index <= curr) {
    stop(sprintf(
      "%s has already run step index %d (to model step %d)",
      name, step_index, steps[step_index, 2]))
  }
}


## This helper pulls together small pieces of bookkeeping that differ
## strongly depending on if the particle filter is a
## nested/multiparameter filter or not.
particle_filter_state_support <- function(nested) {
  if (nested) {
    list(initial = pfs_initial_nested,
         index = pfs_index_nested,
         compare = pfs_compare_nested)
  } else {
    list(initial = pfs_initial_simple,
         index = pfs_index_simple,
         compare = pfs_compare_simple)
  }
}


particle_filter_early_exit <- function(log_likelihood, min_log_likelihood) {
  if (any(log_likelihood == -Inf)) {
    return(TRUE)
  }
  if (length(log_likelihood) == 1) {
    log_likelihood < min_log_likelihood
  } else if (length(min_log_likelihood) == 1) {
    sum(log_likelihood) < min_log_likelihood
  } else {
    all(log_likelihood < min_log_likelihood)
  }
}


## These functions are either "simple" or "nested"
pfs_initial_simple <- function(model, initial, pars, n_particles) {
  state <- initial(model$info(), n_particles, pars)
  if (is.list(state)) {
    stop("Setting 'step' from initial no longer supported")
  }
  state
}


pfs_initial_nested <- function(model, initial, pars, n_particles) {
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


pfs_index_simple <- function(model, index) {
  if (is.null(index)) {
    return(NULL)
  }
  index(model$info())
}


pfs_index_nested <- function(model, index) {
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


pfs_compare_simple <- function(state, compare, data, pars) {
  log_weights <- compare(state, data, pars)
  if (is.null(log_weights)) {
    return(log_weights)
  }
  scale_log_weights(log_weights)
}


pfs_compare_nested <- function(state, compare, data, pars) {
  log_weights <- lapply(seq_len(nlayer(state)), function(i)
    compare(array_drop(state[, , i, drop = FALSE], 3), data[[i]], pars[[i]]))
  if (all(lengths(log_weights) == 0)) {
    return(NULL)
  }

  weights <- lapply(log_weights, scale_log_weights)
  list(average = vnapply(weights, "[[", "average"),
       weights = vapply(weights, function(x) x$weights, numeric(ncol(state))))
}
