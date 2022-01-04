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
    save_restart = NULL
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
      n_particles <- 1L
      if (is.null(model)) {
        model <- generator$new(pars = pars, step = steps[[1]],
                               n_particles = n_particles, n_threads = n_threads,
                               seed = NULL, deterministic = TRUE)
      } else {
        model$update_state(pars = pars, step = steps[[1]])
      }

      if (!is.null(initial)) {
        initial_data <- initial(model$info(), n_particles, pars)
        if (is.list(initial_data)) {
          stop("Setting 'step' from initial no longer supported")
        }
        model$update_state(state = initial_data)
      }

      if (is.null(index)) {
        index_data <- NULL
      } else {
        ## NOTE: this assumes that all parameterisation result in the
        ## same shape, which is assumed generally.
        index_data <- deterministic_index(index(model$info()))
        model$set_index(index_data$index)
      }

      if (save_history) {
        len <- nrow(steps) + 1L
        state <- model$state(index_data$predict)
        history_value <- array(NA_real_, c(dim(state), len))
        history_value[, , 1] <- state
        rownames(history_value) <- names(index_data$predict)
        self$history <- list(
          value = history_value,
          index = index_data$predict)
      }

      save_restart_step <- check_save_restart(save_restart, data)
      if (length(save_restart_step) > 0) {
        self$restart_state <- array(NA_real_,
                                    c(model$n_state(),
                                      n_particles,
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
      private$n_threads <- n_threads
      private$initial <- initial
      private$index <- index
      private$index_data <- index_data
      private$compare <- compare
      private$save_history <- save_history
      private$save_restart_step <- save_restart_step
      private$save_restart <- save_restart

      ## Variable (see also history)
      self$model <- model
      self$log_likelihood <- 0.0
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

      if (save_restart) {
        phases <- deterministic_steps_restart(
          private$save_restart_step, step_end)
        y <- vector("list", length(phases))
        for (i in seq_along(phases)) {
          phase <- phases[[i]]
          y[[i]] <- model$simulate(phase$step_end)
          if (!is.na(phase$restart)) {
            restart_state[, , phase$restart] <- model$state()
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
      } else {
        y_compare <- y[index$run, , , drop = FALSE]
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
        deterministic_likelihood(y_compare, private$compare,
                                 private$pars, private$data_split[idx])

      if (save_history) {
        if (!is.null(index)) {
          y <- y[index$state, , , drop = FALSE]
        }
        ## We need to offset by 1 to allow for the initial conditions
        self$history$value[, , idx + 1L] <- y
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

      info_old <- self$model$info()
      info_new <- ret$model$info()

      state <- transform_state(self$model$state(), info_old, info_new)
      step <- self$model$step()

      ret$model$update_state(state = state, step = step)
      ret$current_step_index <- self$current_step_index
      ret$log_likelihood <- self$log_likelihood

      ret
    }
  ))


deterministic_likelihood <- function(y, compare, pars, data) {
  ll <- numeric(length(data))
  for (i in seq_along(ll)) {
    y_i <- array_drop(y[, 1L, i, drop = FALSE], 3L)
    ll[i] <- compare(y_i, data[[i]], pars)
  }
  sum(ll)
}


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


deterministic_index <- function(index) {
  index_all <- union(index$run, index$state)
  list(index = index_all,
       run = set_names(match(index$run, index_all), names(index$run)),
       state = set_names(match(index$state, index_all), names(index$state)),
       predict = index$state)
}
