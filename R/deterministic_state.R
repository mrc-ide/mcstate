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
    steps = NULL,
    n_threads = NULL,
    initial = NULL,
    index = NULL,
    index_data = NULL,
    compare = NULL,
    save_history = NULL,
    save_restart_step = NULL,
    save_restart = NULL,
    support = NULL
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

    ##' @field current_step_index The index of the last completed step.
    current_step_index = 0L,

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
    ##' @param steps A matrix of step beginning and ends
    ##' @param n_threads The number of threads to use
    ##' @param initial Initial condition function (or `NULL`)
    ##' @param index Index function (or `NULL`)
    ##' @param compare Compare function
    ##' @param save_history Logical, indicating if we should save history
    ##' @param save_restart Vector of steps to save restart at
    initialize = function(pars, generator, model, data, data_split, steps,
                          n_threads, initial, index, compare,
                          save_history, save_restart) {
      pars_multi <- inherits(data, "particle_filter_data_nested")
      support <- particle_deterministic_state_support(pars_multi)

      ## This adds an extra dimension (vs using NULL), which is not
      ## amazing, but it does simplify logic in a few places and keeps
      ## this behaving more similarly to the particle filter.
      n_particles <- 1L
      if (is.null(model)) {
        model <- generator$new(pars = pars, step = steps[[1]],
                               n_particles = n_particles, n_threads = n_threads,
                               seed = NULL, deterministic = TRUE,
                               pars_multi = pars_multi)
      } else {
        model$update_state(pars = pars, step = steps[[1]])
      }

      if (!is.null(initial)) {
        initial_data <- support$initial(model, initial, pars, n_particles)
        model$update_state(state = initial_data)
      }

      if (is.null(index)) {
        index_data <- NULL
      } else {
        index_data <- support$index(model, index)
        model$set_index(index_data$index)
      }

      ## The model shape is [n_particles, <any multi-par structure>]
      shape <- model$shape()

      if (save_history) {
        len <- nrow(steps) + 1L
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
      private$n_threads <- n_threads
      private$initial <- initial
      private$index <- index
      private$index_data <- index_data
      private$compare <- compare
      private$save_history <- save_history
      private$save_restart_step <- save_restart_step
      private$save_restart <- save_restart
      private$support <- support

      ## Variable (see also history)
      self$model <- model
      self$log_likelihood <- rep(0, prod(shape[-1]))
    },

    ##' @description Run the deterministic particle to the end of the data.
    ##' This is a convenience function around `$step()` which provides the
    ##' correct value of `step_index`
    run = function() {
      self$step(nrow(private$steps))
    },

    ##' @description Take a step with the deterministic particle. This moves
    ##' the system forward one step within the *data* (which
    ##' may correspond to more than one step with your model) and
    ##' returns the likelihood so far.
    ##'
    ##' @param step_index The step *index* to move to. This is not the same
    ##' as the model step, nor time, so be careful (it's the index within
    ##' the data provided to the filter). It is an error to provide
    ##' a value here that is lower than the current step index, or past
    ##' the end of the data.
    step = function(step_index) {
      curr <- self$current_step_index
      check_step(curr, step_index, private$steps)

      model <- self$model
      index <- private$index_data

      save_history <- private$save_history

      restart_state <- self$restart_state
      save_restart <- !is.null(restart_state)

      ## Unlike the normal particle filter, we do this all in one shot
      idx <- (curr + 1):step_index
      step_end <- private$steps[idx, 2]

      support <- private$support

      if (save_restart) {
        phases <- deterministic_steps_restart(
          private$save_restart_step, step_end)
        y <- vector("list", length(phases))
        for (i in seq_along(phases)) {
          phase <- phases[[i]]
          y[[i]] <- model$simulate(phase$step_end)
          if (!is.na(phase$restart)) {
            array_last_dimension(restart_state, phase$restart) <- model$state()
          }
        }
        self$restart_state <- restart_state
        y <- array_bind(arrays = y)
      } else {
        y <- model$simulate(step_end)
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
      self$current_step_index <- step_index

      log_likelihood
    },

    ##' @description Create a new `deterministic_particle_state` object based
    ##' on this one (same model, position in time within the data) but with
    ##' new parameters, to support the "multistage particle filter".
    ##'
    ##' @param pars New model parameters
    ##'
    ##' @param transform_state A function to transform the model state
    ##'   from the old to the new parameter set.  See
    ##'   [mcstate::multistage_epoch()] for details.
    fork_multistage = function(pars, transform_state) {
      model <- NULL
      save_history <- !is.null(self$history)
      initial <- NULL

      if (is.null(pars)) {
        pars <- self$model$pars()
      }
      ret <- particle_deterministic_state$new(
        pars, private$generator, model, private$data, private$data_split,
        private$steps, private$n_threads, initial, private$index,
        private$compare, save_history, private$save_restart)

      particle_filter_update_state(transform_state, self$model, ret$model)

      ret$current_step_index <- self$current_step_index
      ret$log_likelihood <- self$log_likelihood

      ret
    }
  ))


deterministic_steps_restart <- function(save_restart_step, step_end) {
  i <- match(save_restart_step, step_end)
  j <- which(!is.na(i))

  ## No restart in this block, do the easy exit:
  if (length(j) == 0L) {
    return(list(list(step_end = step_end, restart = NA_integer_)))
  }

  i <- i[j]
  if (length(i) == 0 || last(i) < length(step_end)) { # first part now dead?
    i <- c(i, length(step_end))
    j <- c(j, NA_integer_)
  }

  ## This feels like it could be done more efficiently this is at
  ## least fairly compact:
  step_end <- unname(split(step_end, rep(seq_along(i), diff(c(0, i)))))
  Map(list, step_end = step_end, restart = j)
}


particle_deterministic_state_support <- function(nested) {
  if (nested) {
    list(initial = pfs_initial_nested,
         index = pds_index_nested,
         compare = pds_compare_nested)
  } else {
    list(initial = pfs_initial_simple,
         index = pds_index_simple,
         compare = pds_compare_simple)
  }
}


pds_index_simple <- function(model, index) {
  pds_index_process(index(model$info()))
}


pds_index_nested <- function(model, index) {
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


pds_compare_simple <- function(state, compare, data, pars) {
  n_data <- length(data)
  ll <- numeric(n_data)
  for (i in seq_len(n_data)) {
    data_i <- data[[i]]
    state_i <- array_drop(state[, 1L, i, drop = FALSE], 3L)
    ll[i] <- compare(state_i, data_i, pars)
  }
  sum(ll)
}


pds_compare_nested <- function(state, compare, data, pars) {
  ## At this point, I really have no strong idea what is going on with
  ## the data!
  n_pars <- length(pars)
  n_data <- length(data)
  ll <- array(0, c(n_data, n_pars))
  for (i in seq_len(n_data)) {
    data_i <- data[[i]]
    state_i <- array_drop(state[, 1L, , i, drop = FALSE], 4)
    for (j in seq_len(n_pars)) {
      ## Drop dimension 3 (parameter) and 4 (time) to leave 1 (state)
      ## and 2 (particle) for compatibility with the particle filter
      ## comparison functions.
      state_ij <- array_drop(state[, 1L, j, i, drop = FALSE], c(3, 4))
      ll[i, j] <- compare(state_ij, data_i[[j]], pars[[j]])
    }
  }
  colSums(ll)
}
