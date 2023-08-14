##' @title Deterministic particle state
##'
##' @description Deterministic particle internal state. This object is
##'   not ordinarily constructed directly by users, but via the
##'   `$run_begin` method to [mcstate::particle_deterministic]. It
##'   provides an advanced interface to the deterministic particle
##'   that allows partially running over part of the time trajectory.
##'
##' This state object has a number of public fields that you can read
##'   but must not write (they are not read-only so you *could* write
##'   them, but don't).
particle_deterministic_state <- R6::R6Class(
  "particle_deterministic_state",
  cloneable = FALSE,

  ## as for particle_filter_state, but missing: n_particles (always 1
  ## per parameter set), gpu (never allowed) and min_log_likelihood
  ## (not useful here)
  private = list(
    generator = NULL,
    pars = NULL,
    data = NULL,
    data_split = NULL,
    times = NULL,
    n_threads = NULL,
    has_multiple_parameters = NULL,
    initial = NULL,
    index = NULL,
    index_data = NULL,
    compare = NULL,
    save_history = NULL,
    save_restart_time = NULL,
    save_restart = NULL,
    support = NULL,

    step_r = function(time_index) {
      curr <- self$current_time_index
      check_time_step(curr, time_index, private$times)

      model <- self$model
      index <- private$index_data

      save_history <- private$save_history

      restart_state <- self$restart_state
      save_restart <- !is.null(restart_state)

      ## Unlike the normal particle filter, we do this all in one shot
      idx <- (curr + 1):time_index
      time_end <- private$times[idx, 2]

      support <- private$support

      if (save_restart) {
        phases <- deterministic_times_restart(
          private$save_restart_time, time_end)
        y <- vector("list", length(phases))
        for (i in seq_along(phases)) {
          phase <- phases[[i]]
          y[[i]] <- model$simulate(phase$time_end)
          if (!is.na(phase$restart)) {
            array_last_dimension(restart_state, phase$restart) <- model$state()
          }
        }
        self$restart_state <- restart_state
        y <- array_bind(arrays = y)
      } else {
        y <- model$simulate(time_end)
        restart_state <- NULL
      }

      if (is.null(index)) {
        y_compare <- y
        y_history <- y
      } else {
        y_compare <- array_first_dimension(y, index$run)
        y_history <- array_first_dimension(y, index$state)
        rownames(y_compare) <- names(index$run)
      }

      ## The likelihood is a loop over both:
      ##
      ## 1. the different parameter sets (because these might affect
      ##    the observation function)
      ## 2. the different timesteps which have different data points
      ##
      ## A clever function might be able to vecorise away one or both
      ## of these.  We can immediately run this over all parameter
      ## sets at once if we do not have parameters involved in the
      ## observation function too.
      log_likelihood <- self$log_likelihood +
        support$compare(y_compare, private$compare, private$data_split[idx],
                        private$pars)

      if (save_history) {
        ## We need to offset by 1 to allow for the initial conditions
        array_last_dimension(self$history$value, idx + 1L) <- y_history
      } else {
        self$history <- NULL
      }

      self$log_likelihood <- log_likelihood
      self$current_time_index <- time_index

      log_likelihood
    },

    step_compiled = function(time_index) {
      curr <- self$current_time_index
      check_time_step(curr, time_index, private$times, "Particle filter")
      time <- private$times[time_index, 2]

      model <- self$model

      history <- self$history
      save_history <- !is.null(history)

      res <- model$filter(time, save_history, private$save_restart_time)

      self$log_likelihood <- self$log_likelihood + res$log_likelihood
      self$current_time_index <- time_index
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
    ##'   `save_history = TRUE`.
    history = NULL,

    ##' @field restart_state Full model state at a series of points in
    ##'   time, if the model was created with non-`NULL` `save_restart`.
    ##'   This is a 3d array as described in [mcstate::particle_filter]
    restart_state = NULL,

    ##' @field log_likelihood The log-likelihood so far. This starts at
    ##'   0 when initialised and accumulates value for each step taken.
    log_likelihood = NULL,

    ##' @field current_time_index The index of the last completed step.
    current_time_index = 0L,

    ## As for private fields; missing
    ## n_particles, gpu_config, min_log_likelihood but also missing seed
    ##' @description Initialise the deterministic particle state. Ordinarily
    ##' this should not be called by users, and so arguments are barely
    ##' documented.
    ##'
    ##' @param pars Parameters for a single phase
    ##' @param generator A dust generator object
    ##' @param model If the generator has previously been initialised
    ##' @param data A [mcstate::particle_filter_data] data object
    ##' @param data_split The same data as `data` but split by step
    ##' @param times A matrix of time step beginning and ends
    ##' @param has_multiple_parameters Compute multiple likelihoods at once?
    ##' @param n_threads The number of threads to use
    ##' @param initial Initial condition function (or `NULL`)
    ##' @param index Index function (or `NULL`)
    ##' @param compare Compare function
    ##' @param constant_log_likelihood Constant log likelihood function
    ##' @param save_history Logical, indicating if we should save history
    ##' @param save_restart Vector of time steps to save restart at
    ##' @param stochastic_schedule Vector of times to perform stochastic updates
    ##' @param ode_control Tuning control for stepper
    initialize = function(pars, generator, model, data, data_split, times,
                          has_multiple_parameters, n_threads,
                          initial, index, compare,
                          constant_log_likelihood,
                          save_history, save_restart,
                          stochastic_schedule, ode_control) {
      has_multiple_data <- inherits(data, "particle_filter_data_nested")
      is_continuous <- inherits(data, "particle_filter_data_continuous")
      support <- particle_deterministic_state_support(has_multiple_parameters,
                                                      has_multiple_data)

      ## This adds an extra dimension (vs using NULL), which is not
      ## amazing, but it does simplify logic in a few places and keeps
      ## this behaving more similarly to the particle filter.
      n_particles <- 1L
      if (is.null(model)) {
        if (is_continuous) {
          model <- generator$new(pars = pars, time = times[[1L]],
                                 n_particles = n_particles,
                                 n_threads = n_threads,
                                 seed = NULL, deterministic = TRUE,
                                 ode_control = ode_control,
                                 pars_multi = has_multiple_parameters)
          model$set_stochastic_schedule(stochastic_schedule)
        } else {
          model <- generator$new(pars = pars, time = times[[1L]],
                                 n_particles = n_particles,
                                 n_threads = n_threads,
                                 seed = NULL, deterministic = TRUE,
                                 pars_multi = has_multiple_parameters)
        }
        if (is.null(compare)) {
          data_is_shared <- has_multiple_parameters && !has_multiple_data
          model$set_data(data_split, data_is_shared)
        }
      } else {
        model$update_state(pars = pars, time = times[[1]])
      }

      if (!is.null(initial)) {
        initial_data <- support$initial(model, initial, pars, n_particles)
        model$update_state(state = initial_data)
      }

      if (is.null(index)) {
        index_data <- NULL
      } else {
        index_data <- support$index(model, index)
        if (is.null(compare)) {
          model$set_index(index_data$predict)
        } else {
          model$set_index(index_data$index)
        }
      }

      ## The model shape is [n_particles, <any multi-par structure>]
      shape <- model$shape()

      if (save_history) {
        len <- nrow(times) + 1L
        state <- model$state(index_data$predict)
        history_value <- array(NA_real_, c(dim(state), len))
        array_last_dimension(history_value, 1) <- state
        rownames(history_value) <- names(index_data$predict)
        self$history <- list(
          value = history_value,
          index = index_data$predict)
      } else {
        self$history <- NULL
      }

      save_restart_time <- check_save_restart(save_restart, data)
      if (length(save_restart_time) > 0) {
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
      private$has_multiple_parameters <- has_multiple_parameters
      private$n_threads <- n_threads
      private$initial <- initial
      private$index <- index
      private$index_data <- index_data
      private$compare <- compare
      private$save_history <- save_history
      private$save_restart_time <- save_restart_time
      private$save_restart <- save_restart
      private$support <- support

      ## Variable (see also history)
      self$model <- model
      self$log_likelihood <- particle_filter_constant_log_likelihood(
        pars, has_multiple_parameters, constant_log_likelihood)
    },

    ##' @description Run the deterministic particle to the end of the data.
    ##' This is a convenience function around `$step()` which provides the
    ##' correct value of `time_index`
    run = function() {
      self$step(nrow(private$times))
    },

    ##' @description Take a step with the deterministic particle. This moves
    ##' the system forward one step within the *data* (which
    ##' may correspond to more than one step with your model) and
    ##' returns the likelihood so far.
    ##'
    ##' @param time_index The step *index* to move to. This is not the same
    ##' as the model step, nor time, so be careful (it's the index within
    ##' the data provided to the filter). It is an error to provide
    ##' a value here that is lower than the current step index, or past
    ##' the end of the data.
    step = function(time_index) {
      if (is.null(private$compare)) {
        private$step_compiled(time_index)
      } else {
        private$step_r(time_index)
      }
    },

    ##' @description Create a new `deterministic_particle_state` object based
    ##' on this one (same model, position in time within the data) but with
    ##' new parameters, to support the "multistage particle filter".
    ##'
    ##' @param model A model object
    ##'
    ##' @param pars New model parameters
    ##'
    ##' @param transform_state A function to transform the model state
    ##'   from the old to the new parameter set.  See
    ##'   [mcstate::multistage_epoch()] for details.
    fork_multistage = function(model, pars, transform_state) {
      save_history <- !is.null(self$history)
      initial <- NULL
      constant_log_likelihood <- NULL

      if (is.null(pars)) {
        pars <- self$model$pars()
      }
      ret <- particle_deterministic_state$new(
        pars, private$generator, model, private$data, private$data_split,
        private$times, private$has_multiple_parameters, private$n_threads,
        initial, private$index, private$compare, constant_log_likelihood,
        save_history, private$save_restart)

      particle_filter_update_state(transform_state, self$model, ret$model)

      ret$current_time_index <- self$current_time_index
      ret$log_likelihood <- self$log_likelihood

      ret
    }
  ))


deterministic_times_restart <- function(save_restart_time, time_end) {
  i <- match(save_restart_time, time_end)
  j <- which(!is.na(i))

  ## No restart in this block, do the easy exit:
  if (length(j) == 0L) {
    return(list(list(time_end = time_end, restart = NA_integer_)))
  }

  i <- i[j]
  if (length(i) == 0 || last(i) < length(time_end)) { # first part now dead?
    i <- c(i, length(time_end))
    j <- c(j, NA_integer_)
  }

  ## This feels like it could be done more efficiently this is at
  ## least fairly compact:
  time_end <- unname(split(time_end, rep(seq_along(i), diff(c(0, i)))))
  Map(list, time_end = time_end, restart = j)
}


particle_deterministic_state_support <- function(has_multiple_parameters,
                                                 has_multiple_data) {
  if (has_multiple_parameters) {
    list(initial = pfs_initial_multiple,
         index = pds_index_multiple,
         compare = function(...) pds_compare_multiple(has_multiple_data, ...))
  } else {
    list(initial = pfs_initial_single,
         index = pds_index_single,
         compare = pds_compare_single)
  }
}


pds_index_single <- function(model, index) {
  pds_index_process(index(model$info()))
}


pds_index_multiple <- function(model, index) {
  index_data <- lapply(model$info(), index)

  nok <- !all(vlapply(index_data[-1], identical, index_data[[1]]))
  if (nok) {
    stop("index must be identical across populations")
  }

  pds_index_process(index_data[[1]])
}


pds_index_process <- function(data) {
  index_all <- union(data$run, data$state)
  list(index = index_all,
       run = set_names(match(data$run, index_all), names(data$run)),
       state = set_names(match(data$state, index_all), names(data$state)),
       predict = data$state)
}


pds_compare_single <- function(state, compare, data, pars) {
  n_data <- length(data)
  ll <- numeric(n_data)
  for (i in seq_len(n_data)) {
    data_i <- data[[i]]
    state_i <- array_drop(state[, 1L, i, drop = FALSE], 3L)
    ll[i] <- compare(state_i, data_i, pars)
  }
  sum(ll)
}


pds_compare_multiple <- function(has_multiple_data, state, compare, data,
                                 pars) {
  n_pars <- length(pars) # number of parameter sets
  n_data <- length(data) # number of time points in data (*not* data sets)
  ll <- array(0, c(n_data, n_pars))
  for (i in seq_len(n_data)) {
    data_i <- data[[i]]
    state_i <- array_drop(state[, 1L, , i, drop = FALSE], 4)
    for (j in seq_len(n_pars)) {
      ## Drop dimension 3 (parameter) and 4 (time) to leave 1 (state)
      ## and 2 (particle) for compatibility with the particle filter
      ## comparison functions.
      state_ij <- array_drop(state[, 1L, j, i, drop = FALSE], c(3, 4))
      data_ij <- if (has_multiple_data) data_i[[j]] else data_i
      ll[i, j] <- compare(state_ij, data_ij, pars[[j]])
    }
  }
  colSums(ll)
}
