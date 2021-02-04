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
    save_restart_step = NULL,
    save_restart = NULL,
    current_step_index = 0L
  ),

  public = list(
    ##' @field model The dust model generator being simulated (cannot be
    ##'   re-bound)
    model = NULL,

    ##' @field complete Logical, indicating if the particle filter has
    ##'   reached the end of the data.
    complete = FALSE,

    ##' @field history The particle history, if created with
    ##'   `save_history = TRUE`. This is an internal format subject to
    ##    change.
    history = NULL,

    ##' @field restart_state Full model state at a series of points in
    ##'   time, if the model was created with non-`NULL` `save_restart`.
    ##'   This is an internal format subject to change.
    restart_state = NULL,

    ##' @field log_likelihood The log-likelihood so far. This starts at
    ##'   0 when initialised and accumulates value for each step taken.
    log_likelihood = NULL,

    ##' @field index_state The index used to query state when running
    ##'   the model (this is needed to make sense of the final model state.
    ## TODO: I think this can be removed.
    index_state = NULL,

    ##' @description Initialise the particle filter state. Ordinarily
    ##' this should not be called by users, and so arguments are not
    ##' documented.
    initialize = function(pars, generator, model, data, data_split, steps,
                          n_particles, n_threads, initial, index, compare,
                          seed, save_history, save_restart) {
      ## NOTE: this will generate a warning when updating docs but
      ## that's ok; see https://github.com/r-lib/roxygen2/issues/1067
      if (is.null(model)) {
        model <- generator$new(pars = pars, step = steps[[1L]],
                               n_particles = n_particles, n_threads = n_threads,
                               seed = seed)
        if (is.null(compare)) {
          model$set_index(integer(0))
          model$set_data(data_split)
        }
      } else {
        model$reset(pars, steps[[1L]])
      }

      if (!is.null(initial)) {
        initial_data <- initial(model$info(), n_particles, pars)
        if (is.list(initial_data)) {
          steps <- particle_steps(steps, initial_data$step)
          model$set_state(initial_data$state, initial_data$step)
        } else {
          model$set_state(initial_data)
        }
      }

      index_data <- if (is.null(index)) NULL else index(model$info())
      if (!is.null(compare) && !is.null(index_data$run)) {
        model$set_index(index_data$run)
      }

      if (save_history) {
        len <- nrow(steps) + 1L
        state <- model$state(index_data$state)
        history_value <- array(NA_real_, c(dim(state), len))
        history_value[, , 1] <- state
        history_order <- matrix(seq_len(n_particles), n_particles, len)
        self$history <- list(
          value = history_value,
          order = history_order,
          index = index_data$state)
        ## We also need to export this out because the *final state*
        ## needs to know this.
        self$index_state <- index_data$state
      } else {
        self$history <- NULL
      }

      ## TODO: replace this with an array, once we know all the
      ## dimensions here; that requires exposing n_state from the dust
      ## object to do efficiently and generally.
      save_restart_step <- check_save_restart(save_restart, data)
      if (length(save_restart_step) > 0) {
        self$restart_state <- history_collector(length(save_restart_step) - 1L)
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
      private$save_restart_step <- save_restart_step
      private$save_restart <- save_restart

      ## Variable (see also history)
      self$model <- model
      self$log_likelihood <- 0.0
    },

    ##' @description Take a step with the particle filter. This moves
    ##' the particle filter forward one step within the *data* (which
    ##' may correspond to more than one step with your model) and
    ##' returns the likelihood so far.
    step = function() {
      if (self$complete) {
        stop("The particle filter has reached the end of the data")
      }
      private$current_step_index <- private$current_step_index + 1L
      self$complete <- private$current_step_index >= nrow(private$steps)
      step <- private$current_step_index
      step_end <- private$steps[step, 2]
      save_history <- !is.null(self$history)

      state <- self$model$run(step_end)
      if (save_history) {
        self$history$value[, , step + 1L] <-
          self$model$state(self$history$index)
      }

      if (is.null(private$compare)) {
        log_weights <- self$model$compare_data()
      } else {
        log_weights <- private$compare(state, private$data_split[[step]],
                                       private$pars)
      }

      if (is.null(log_weights)) {
        if (save_history) {
          self$history$order[, step + 1L] <- seq_len(ncol(state))
        }
      } else {
        weights <- scale_log_weights(log_weights)
        self$log_likelihood <- self$log_likelihood + weights$average
        if (self$log_likelihood == -Inf) {
          self$complete <- TRUE
          return(-Inf)
        }

        kappa <- particle_resample(weights$weights)
        self$model$reorder(kappa)
        if (save_history) {
          self$history$order[, step + 1L] <- kappa
        }
      }

      if (step_end %in% private$save_restart_step) {
        self$restart_state$add(self$model$state())
      }

      self$log_likelihood
    },

    ##' @description Create a new `particle_filter_state` object based
    ##' on this one (same model, position in time within the data) but
    ##' with new parameters. To do this, we create a new
    ##' `particle_filter_state` with new parameters at the beginning of
    ##' the simulation (corresponding to the start of your data or the
    ##' `initial` argument to [mcstate::particle_filter]) with your new
    ##' `pars`, and then run the filter foward in time until it reaches
    ##' the same step as the parent model.
    ##'
    ##' @param pars New model parameters
    fork = function(pars) {
      seed <- self$model$rng_state()
      save_history <- !is.null(self$history)
      ret <- particle_filter_state$new(
        pars, private$generator, NULL, private$data, private$data_split,
        private$steps, private$n_particles, private$n_threads,
        private$initial, private$index, private$compare, seed,
        save_history, private$save_restart)

      ## Run it up to the same point
      for (i in seq_len(private$current_step_index)) {
        ret$step()
      }

      ## Set the seed in the parent model
      self$model$set_rng_state(ret$model$rng_state())

      ## Now, the user can use this model
      ret
    }
  ))
