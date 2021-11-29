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
    n_steps = NULL,
    n_particles = NULL,
    n_threads = NULL,
    initial = NULL,
    index = NULL,
    compare = NULL,
    gpu = NULL,
    save_restart_step = NULL,
    save_restart = NULL,
    min_log_likelihood = NULL,

    spawn = function(pars, state) {
      stopifnot(!private$gpu) # this won't work
      gpu_config <- NULL
      seed <- self$model$rng_state()
      save_history <- !is.null(self$history)
      if (is.null(state)) {
        initial <- private$initial
      } else {
        initial <- function(...) {
          list(state = state, step = self$model$step())
        }
      }
      model <- NULL
      particle_filter_state$new(
        pars, private$generator, model, private$data, private$data_split,
        private$steps, private$n_particles, private$n_threads,
        initial, private$index, private$compare, gpu_config,
        seed, private$min_log_likelihood, save_history, private$save_restart)
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
    ##'   This is a 3d array as described in [mcstate::particle_filter]
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
    ##' this should not be called by users, and so arguments are not
    ##' documented.
    initialize = function(pars, generator, model, data, data_split, steps,
                          n_particles, n_threads, initial, index, compare,
                          gpu_config, seed, min_log_likelihood,
                          save_history, save_restart) {
      ## NOTE: this will generate a warning when updating docs but
      ## that's ok; see https://github.com/r-lib/roxygen2/issues/1067
      if (is.null(model)) {
        model <- generator$new(pars = pars, step = steps[[1L]],
                               n_particles = n_particles, n_threads = n_threads,
                               seed = seed, gpu_config = gpu_config)
        if (is.null(compare)) {
          model$set_index(integer(0))
          model$set_data(data_split)
        }
      } else {
        model$update_state(pars = pars, step = steps[[1L]])
      }

      if (!is.null(initial)) {
        initial_data <- initial(model$info(), n_particles, pars)
        if (is.list(initial_data)) {
          steps <- particle_steps(steps, initial_data$step)
          model$update_state(state = initial_data$state,
                             step = initial_data$step)
        } else {
          model$update_state(state = initial_data)
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
      } else {
        self$history <- NULL
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
      private$n_steps <- nrow(steps)
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
      self$log_likelihood <- 0.0
    },

    ##' @description Run the particle filter to the end of the data. This is
    ##' a convenience function around `$step()` which provides the correct
    ##' value of `step_index`
    run = function() {
      if (is.null(private$compare)) {
        ## TODO: add min_log_likelihood support here, needs work in dust?
        particle_filter_compiled(self, private)
      } else {
        self$step(private$n_steps)
      }
    },

    continue = function(other, state) {
      self$model$set_rng_state(other$model$rng_state())
      self$update_state(other$current_step_index, state)
      self$log_likelihood_step <- other$log_likelihood_step
      self$log_likelihood <- other$log_likelihood
    },

    update_state = function(step_index, state) {
      self$current_step_index <- step_index
      step <- private$steps[step_index, 2]
      self$model$update_state(state = state, step = step)
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
      n_steps <- private$n_steps
      curr <- self$current_step_index
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
          history_value[, , t + 1L] <- model$state(save_history_index)
        }

        log_weights <- compare(state, data_split[[t]], pars)

        if (is.null(log_weights)) {
          if (save_history) {
            history_order[, t + 1L] <- seq_len(n_particles)
          }
          log_likelihood_step <- NA_real_
        } else {
          weights <- scale_log_weights(log_weights)
          log_likelihood_step <- weights$average
          log_likelihood <- log_likelihood + log_likelihood_step
          if (log_likelihood < min_log_likelihood) {
            log_likelihood <- -Inf
          }
          if (log_likelihood == -Inf) {
            break
          }

          kappa <- particle_resample(weights$weights)
          model$reorder(kappa)
          if (save_history) {
            history_order[, t + 1L] <- kappa
          }
        }

        if (save_restart) {
          i_restart <- match(step_end, save_restart_step)
          if (!is.na(i_restart)) {
            restart_state[, , i_restart] <- model$state()
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

    update = function(pars, state) {
      ## this is *very* similar to fork, which exists to support pmc2, except:
      ## * we provide state
      ## * we could continue with the same model object, though this is
      ##   not implemented yet.
      private$spawn(pars, state)
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
      ret <- private$spawn(pars, NULL)

      ## Run it up to the same point
      ret$step(self$current_step_index)

      ## Set the seed in the parent model
      self$model$set_rng_state(ret$model$rng_state())

      ## Now, the user can use this model
      ret
    }
  ))


## This is used by both the nested and non-nested particle filter, and
## outsources all the work to dust.
particle_filter_compiled <- function(self, private) {
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
  self$current_step_index <- private$n_steps
  if (save_history) {
    self$history <- list(value = res$trajectories,
                         index = save_history_index)
  }
  self$restart_state <- res$snapshots
  res$log_likelihood
}
