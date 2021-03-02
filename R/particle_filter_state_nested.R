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
    save_restart_step = NULL,
    save_restart = NULL,
    current_step_index = 0L
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

    ##' @description Initialise the particle filter state. Ordinarily
    ##' this should not be called by users, and so arguments are not
    ##' documented.
    initialize = function(pars, generator, model, data, data_split, steps,
                          n_particles, n_threads, initial, index, compare,
                          seed, save_history, save_restart) {

      ## NOTE: this will generate a warning when updating docs but
      ## that's ok; see https://github.com/r-lib/roxygen2/issues/1067
      if (length(pars) < length(unique(data$population))) {
        stop(sprintf("'pars' must be at least the length of
                     'data$population', %d",
                     length(unique(data$population))))
      }

      if (is.null(model)) {
        model <- generator$new(pars = pars, step = steps[[1L]],
                               n_particles = n_particles,
                               n_threads = n_threads, seed = seed,
                               pars_multi = TRUE)
        if (is.null(compare)) {
          model$set_index(integer(0))
          model$set_data(data_split)
        }
      } else {
        model$reset(pars, steps[[1L]])
      }

      if (!is.null(initial)) {
        initial_data <- Map(initial,
          info = model$info(),
          pars = pars,
          MoreArgs = list(n_particles = n_particles)
        )
        if (!is.null(unlist(initial_data))) {
          if (any(c(c("step", "state") %in% names(initial_data[[1]])))) {
            init_step <- unlist(lapply(initial_data, "[[", "step"))
            if (length(init_step) != n_particles * length(pars)) {
              if (length(init_step) == length(pars)) {
                init_step <- rep(init_step, each = n_particles)
              } else {
                stop(sprintf("Expected 'step' to be length %d or %d",
                             length(pars), n_particles * length(pars)))
              }
            }
            init_state <- unlist(lapply(initial_data, "[[", "state"))
            if (!is.null(init_state)) {
              init_state <- matrix(init_state, ncol = 2)
            }
            steps <- particle_steps(steps, init_step)
            model$set_state(init_state, init_step)
          } else {
            model$set_state(list_to_array(initial_data))
          }
        } else {
          model$set_state(NULL)
        }
      }

      index_data <- if (is.null(index)) NULL else lapply(model$info(), index)
      if (!is.null(compare) && !is.null(index_data[[1]]$run)) {
        run <- index_data[[1]]$run
        ok <- vlapply(index_data[-1], function(x) identical(x$run, run))
        if (!all(ok)) {
          stop("Populations should have the same run indices")
        }
        model$set_index(run)
      }

      if (save_history) {
        len <- nrow(steps) + 1L
        index_state <- index_data[[1]]$state
        state <- model$state(index_state)
        history_value <- array(NA_real_, c(dim(state), len))
        history_value[, , , 1] <- state
        history_order <- array(seq_len(n_particles),
                               c(n_particles, nlayer(state), len))

        self$history <- list(
          value = history_value,
          order = history_order,
          index = index_state)
      } else {
        self$history <- NULL
      }

      save_restart_step <- check_save_restart(save_restart, data)
      if (length(save_restart_step) > 0) {
        self$restart_state <- array(NA_real_,
                                    c(model$n_state(),
                                      n_particles,
                                      length(pars),
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
      private$n_steps <- nrow(steps)
      private$n_particles <- n_particles
      private$n_threads <- n_threads
      private$initial <- initial
      private$index <- index
      private$compare <- compare
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
      self$step(private$n_steps)
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
      n_steps <- private$n_steps
      curr <- private$current_step_index
      if (curr >= n_steps) {
        stop("Particle filter has reached the end of the data")
      }
      if (step_index > n_steps) {
        stop(sprintf("step_index %d is beyond the length of the data (max %d)",
                     step_index, n_steps))
      }
      if (step_index <= curr) {
        stop(sprintf(
          "Particle filter has already run step index %d (to model step %d)",
          step_index, steps[step_index, 2]))
      }

      model <- self$model

      compare <- private$compare
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

      for (t in seq(curr + 1L, step_index)) {
        step_end <- steps[t, 2L]
        state <- model$run(step_end)

        if (save_history) {
          history_value[, , , t + 1L] <- model$state(save_history_index)
        }

        if (is.null(compare)) {
          log_weights <- model$compare_data()
        } else {
          ok <- !is.null(compare(
            array_drop(state[, , 1, drop = FALSE], 3),
            data_split[[t]][[1]], pars[[1]]))
          
          if (!ok) {
            log_weights <- NULL
          } else {
            log_weights <- vapply(
              seq_len(nlayer(state)),
              function(i) {
                compare(array_drop(state[, , i, drop = FALSE], 3),
                        data_split[[t]][[i]], pars[[i]])
              }, numeric(ncol(state)))
          }

        }

        if (is.null(log_weights)) {
          if (save_history) {
            history_order[, , t + 1L] <- seq_len(n_particles)
          }
          log_likelihood_step <- NA_real_
        } else {
          weights <- apply(log_weights, 2, scale_log_weights)
          log_likelihood_step <- vnapply(weights, "[[", "average")
          log_likelihood <- log_likelihood + log_likelihood_step
          ## FIXME - IS THIS CORRECT? ALL VS ANY
          if (any(log_likelihood == -Inf)) {
            break
          }

          kappa <- vapply(weights, function(x) particle_resample(x$weights),
                          numeric(length(weights[[1]]$weights)))
          model$reorder(kappa)

          if (save_history) {
            history_order[, , t + 1L] <- kappa
          }
        }

        if (save_restart) {
          i_restart <- match(step_end, save_restart_step)
          if (!is.na(i_restart)) {
            restart_state[, , , i_restart] <- model$state()
          }
        }
      }

      self$log_likelihood_step <- log_likelihood_step
      self$log_likelihood <- log_likelihood
      private$current_step_index <- step_index
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
    fork = function(pars) {
      seed <- self$model$rng_state()
      save_history <- !is.null(self$history)
      ret <- particle_filter_state_nested$new(
        pars, private$generator, NULL, private$data, private$data_split,
        private$steps, private$n_particles, private$n_threads,
        private$initial, private$index, private$compare, seed,
        save_history, private$save_restart)

      ## Run it up to the same point
      ret$step(private$current_step_index)

      ## Set the seed in the parent model
      self$model$set_rng_state(ret$model$rng_state())

      ## Now, the user can use this model
      ret
    }
  ))
