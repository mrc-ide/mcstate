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
    times = NULL,
    n_particles = NULL,
    n_threads = NULL,
    has_multiple_parameters = NULL,
    initial = NULL,
    index = NULL,
    compare = NULL,
    gpu = NULL,
    save_restart_step = NULL,
    save_restart = NULL,
    min_log_likelihood = NULL,
    support = NULL,

    step_r = function(step_index, partial = FALSE) {
      times <- private$times
      curr <- self$current_step_index
      check_step(curr, step_index, private$times, "Particle filter")

      model <- self$model
      compare <- private$compare
      has_multiple_parameters <- private$has_multiple_parameters

      times <- private$times
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
        time_end <- times[t, 2L]
        state <- model$run(time_end)

        if (save_history) {
          ## NOTE: There are two places here (and for the order below)
          ## where we assign trajectories, and we have to do this
          ## without using `array_last_dimension<-`, otherwise there's
          ## a big performance regression due to excessive GC (we make
          ## a copy into the function call, I suspect?)
          ##
          ## An alternative approach here would be to store the
          ## history as a flat array (or access it as such), then
          ## compute what the index would be.  If we know the product
          ## of the first dimensions, this is pretty easy really.
          if (has_multiple_parameters) {
            history_value[, , , t + 1L] <- model$state(save_history_index)
          } else {
            history_value[, , t + 1L] <- model$state(save_history_index)
          }
        }

        weights <- support$compare(state, compare, data_split[[t]], pars)

        if (is.null(weights)) {
          log_likelihood_step <- NA_real_
          kappa <- seq_len(n_particles)
        } else {
          log_likelihood_step <- weights$average
          log_likelihood <- log_likelihood + log_likelihood_step

          if (particle_filter_early_exit(log_likelihood, min_log_likelihood)) {
            log_likelihood[] <- -Inf
            break
          }

          kappa <- particle_resample(weights$weights)
          model$reorder(kappa)
        }

        if (save_history) {
          if (has_multiple_parameters) {
            history_order[, , t + 1L] <- kappa
          } else {
            history_order[, t + 1L] <- kappa
          }
        }

        if (save_restart) {
          i_restart <- match(time_end, save_restart_step)
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

    step_compiled = function(step_index, partial = FALSE) {
      if (partial) {
        stop("'partial' not supported with compiled compare")
      }
      curr <- self$current_step_index
      check_step(curr, step_index, private$times, "Particle filter")
      step <- private$times[step_index, 2]

      model <- self$model

      history <- self$history
      save_history <- !is.null(history)

      res <- model$filter(step, save_history, private$save_restart_step)

      self$log_likelihood_step <- NA_real_
      self$log_likelihood <- self$log_likelihood + res$log_likelihood
      self$current_step_index <- step_index
      if (save_history) {
        self$history <- list(value = res$trajectories,
                             index = self$history$index)
      }
      self$restart_state <- res$snapshots
      self$log_likelihood
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
    ##' @param times A matrix of step beginning and ends
    ##' @param n_particles Number of particles to use
    ##' @param has_multiple_parameters Compute multiple likelihoods at once?
    ##' @param n_threads The number of threads to use
    ##' @param initial Initial condition function (or `NULL`)
    ##' @param index Index function (or `NULL`)
    ##' @param compare Compare function
    ##' @param constant_log_likelihood Constant log likelihood function
    ##' @param gpu_config GPU configuration, passed to `generator`
    ##' @param seed Initial RNG seed
    ##' @param min_log_likelihood Early termination control
    ##' @param save_history Logical, indicating if we should save history
    ##' @param save_restart Vector of times to save restart at
    ##' @param stochastic_schedule Vector of times to perform stochastic updates
    ##' @param ode_control Tuning control for stepper
    initialize = function(pars, generator, model, data, data_split, times,
                          n_particles, has_multiple_parameters,
                          n_threads, initial, index, compare,
                          constant_log_likelihood, gpu_config, seed,
                          min_log_likelihood, save_history, save_restart,
                          stochastic_schedule, ode_control) {
      has_multiple_data <- inherits(data, "particle_filter_data_nested")
      is_continuous <- inherits(data, "particle_filter_data_continuous")

      support <- particle_filter_state_support(has_multiple_parameters,
                                               has_multiple_data)

      if (is.null(model)) {
        if (is_continuous) {
          model <- generator$new(pars = pars, time = times[[1L]],
                                 n_particles = n_particles,
                                 n_threads = n_threads,
                                 seed = seed,
                                 control = ode_control)
          model$set_stochastic_schedule(stochastic_schedule)
        } else {
          model <- generator$new(pars = pars, step = times[[1L]],
                                 n_particles = n_particles,
                                 n_threads = n_threads,
                                 seed = seed, gpu_config = gpu_config,
                                 pars_multi = has_multiple_parameters)
        }
        if (is.null(compare)) {
          data_is_shared <- has_multiple_parameters && !has_multiple_data
          model$set_data(data_split, data_is_shared)
        }
      } else {
        if (is_continuous) {
          model$update_state(pars = pars, time = times[[1L]])
        } else {
          model$update_state(pars = pars, step = times[[1L]])
        }
      }

      if (!is.null(initial)) {
        initial_data <- support$initial(model, initial, pars, n_particles)
        model$update_state(state = initial_data)
      }

      if (is.null(index)) {
        index_data <- NULL
      } else {
        index_data <- support$index(model, index)
        if (!is.null(compare)) {
          model$set_index(index_data$run)
        } else {
          model$set_index(index_data$state)
        }
      }

      ## The model shape is [n_particles, <any multi-par structure>]
      if (is_continuous) {
        shape <- model$n_particles()
      } else {
        shape <- model$shape()
      }

      if (save_history) {
        len <- nrow(times) + 1L
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
        stopifnot(!is_continuous)
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
      private$times <- times
      private$n_particles <- n_particles
      private$n_threads <- n_threads
      private$has_multiple_parameters <- has_multiple_parameters
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
      self$log_likelihood <- particle_filter_constant_log_likelihood(
        pars, has_multiple_parameters, constant_log_likelihood)
    },

    ##' @description Run the particle filter to the end of the data. This is
    ##' a convenience function around `$step()` which provides the correct
    ##' value of `step_index`
    run = function() {
      self$step(nrow(private$times))
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
      if (is.null(private$compare)) {
        private$step_compiled(step_index, partial)
      } else {
        private$step_r(step_index, partial)
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
    ##' @param model A model object (or NULL)
    ##'
    ##' @param pars New model parameters
    ##'
    ##' @param transform_state A function to transform the model state
    ##'   from the old to the new parameter set.  See
    ##'   [mcstate::multistage_epoch()] for details.
    fork_multistage = function(model, pars, transform_state) {
      stopifnot(!private$gpu) # this won't work
      gpu_config <- NULL
      seed <- self$model$rng_state()
      save_history <- !is.null(self$history)
      initial <- NULL
      constant_log_likelihood <- NULL

      if (!is.null(model)) {
        model$set_rng_state(seed)
      }

      if (is.null(pars)) {
        pars <- self$model$pars()
      }
      ret <- particle_filter_state$new(
        pars, private$generator, model, private$data, private$data_split,
        private$times, private$n_particles, private$has_multiple_parameters,
        private$n_threads, initial, private$index, private$compare,
        constant_log_likelihood, gpu_config,
        seed, private$min_log_likelihood, save_history, private$save_restart)

      particle_filter_update_state(transform_state, self$model, ret$model)

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
        private$times, private$n_particles, private$has_multiple_parameters,
        private$n_threads, private$initial, private$index, private$compare,
        private$constant_log_likelihood, gpu_config,
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
check_step <- function(curr, step_index, times, name) {
  n_times <- nrow(times)
  if (curr >= n_times) {
    stop(sprintf("%s has reached the end of the data", name))
  }
  if (step_index > n_times) {
    stop(sprintf("step_index %d is beyond the length of the data (max %d)",
                 step_index, n_times))
  }
  if (step_index <= curr) {
    stop(sprintf(
      "%s has already run step index %d (to model step %d)",
      name, step_index, times[step_index, 2]))
  }
}


## This helper pulls together small pieces of bookkeeping that differ
## strongly depending on if the particle filter is a
## single-/multi-parameter filter or not.
particle_filter_state_support <- function(has_multiple_parameters,
                                          has_multiple_data) {
  if (has_multiple_parameters) {
    list(initial = pfs_initial_multiple,
         index = pfs_index_multiple,
         compare = function(...) pfs_compare_multiple(has_multiple_data, ...))
  } else {
    list(initial = pfs_initial_single,
         index = pfs_index_single,
         compare = pfs_compare_single)
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


particle_filter_constant_log_likelihood <- function(pars,
                                                    has_multiple_parameters,
                                                    constant_log_likelihood) {
  if (is.null(constant_log_likelihood)) {
    if (has_multiple_parameters) {
      rep(0, length(pars))
    } else {
      0
    }
  } else {
    if (has_multiple_parameters) {
      vnapply(pars, constant_log_likelihood)
    } else {
      constant_log_likelihood(pars)
    }
  }
}


particle_filter_update_state <- function(transform, model_old, model_new) {
  info_old <- model_old$info()
  info_new <- model_new$info()
  n_pars <- model_old$n_pars()

  if (n_pars == 0) {
    state <- transform(model_old$state(), info_old, info_new)
  } else {
    state_old <- model_old$state()
    state_new <- lapply(seq_len(n_pars), function(i)
      transform(array_drop(state_old[, , i, drop = FALSE], 3),
                info_old[[i]], info_new[[i]]))
    state <- vapply(state_new, identity, state_new[[1]])
  }

  step <- model_old$step()
  model_new$update_state(state = state, step = step)
}


## These functions are either "single" or "multiple" depending on the
## number of parameters a particle filter is operating with.
pfs_initial_single <- function(model, initial, pars, n_particles) {
  state <- initial(model$info(), n_particles, pars)
  if (is.list(state)) {
    stop("Setting 'step' from initial no longer supported")
  }
  state
}


pfs_initial_multiple <- function(model, initial, pars, n_particles) {
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

  array_from_list(state)
}


pfs_index_single <- function(model, index) {
  index(model$info())
}


pfs_index_multiple <- function(model, index) {
  index_data <- lapply(model$info(), index)

  nok <- !all(vlapply(index_data[-1], identical, index_data[[1]]))
  if (nok) {
    stop("index must be identical across populations")
  }

  index_data[[1]]
}


pfs_compare_single <- function(state, compare, data, pars) {
  log_weights <- compare(state, data, pars)
  if (is.null(log_weights)) {
    return(log_weights)
  }
  scale_log_weights(log_weights)
}


## NOTE: funny argument ordering as we'll partially apply the first
## argument later.
pfs_compare_multiple <- function(has_multiple_data, state, compare, data,
                                 pars) {
  if (has_multiple_data) {
    log_weights <- lapply(seq_len(nlayer(state)), function(i)
      compare(array_drop(state[, , i, drop = FALSE], 3), data[[i]], pars[[i]]))
  } else {
    log_weights <- lapply(seq_len(nlayer(state)), function(i)
      compare(array_drop(state[, , i, drop = FALSE], 3), data, pars[[i]]))
  }
  if (all(lengths(log_weights) == 0)) {
    return(NULL)
  }

  weights <- lapply(log_weights, scale_log_weights)
  list(average = vnapply(weights, "[[", "average"),
       weights = vapply(weights, function(x) x$weights, numeric(ncol(state))))
}
