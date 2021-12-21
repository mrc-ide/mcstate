## It is quite possible that we can merge this with
## particle_filter_state but I am also not sure that will be helpful.
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
    n_steps = NULL,
    n_threads = NULL,
    initial = NULL,
    index = NULL,
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
    ##' @description Initialise the particle filter state. Ordinarily
    ##' this should not be called by users, and so arguments are not
    ##' documented.
    initialize = function(pars, generator, model, data, data_split, steps,
                          n_threads, initial, index, compare,
                          save_history, save_restart) {
      n_particles <- length(pars)
      resize <- !is.null(model) && n_particles != model$n_particles()
      if (is.null(model) || resize) {
        model <- generator$new(pars = pars, step = steps[[1]],
                               n_particles = NULL,
                               n_threads = n_threads,
                               seed = NULL,
                               deterministic = TRUE,
                               pars_multi = TRUE)
      } else {
        model$update_state(pars = pars, step = steps[[1]])
      }

      if (!is.null(initial)) {
        initial_data <- Map(function(p, i) initial(i, 1L, p),
                            pars, model$info())
        if (any(vlapply(initial_data, is.list))) {
          stop("Setting 'step' from initial no longer supported")
        }
        state <- vapply(initial_data, identity, initial_data[[1]])
        model$update_state(state = state)
      }

      if (is.null(index)) {
        index <- NULL
      } else {
        ## NOTE: this assumes that all parameterisation result in the
        ## same shape, which is assumed generally.
        index <- deterministic_index(index(model$info()[[1L]]))
        model$set_index(index$index)
      }

      ## Constants
      private$generator <- generator
      private$pars <- pars
      private$data <- data
      private$data_split <- data_split
      private$steps <- steps
      private$n_steps <- nrow(steps)
      private$n_threads <- n_threads
      private$initial <- initial
      private$index <- index
      private$compare <- compare
      private$save_history <- save_history
      ## private$save_restart_step <- save_restart_step
      private$save_restart <- save_restart

      ## Variable (see also history)
      self$model <- model
      self$log_likelihood <- 0.0
    },

    ##' @description Run the particle filter to the end of the data. This is
    ##' a convenience function around `$step()` which provides the correct
    ##' value of `step_index`
    run = function() {
      self$step(private$n_steps)
    },

    step = function(step_index) {
      steps <- private$steps
      n_steps <- private$n_steps
      curr <- self$current_step_index
      ## TODO: Share with particle filter state
      if (curr >= n_steps) {
        stop("Particle deterministic has reached the end of the data")
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
      index <- private$index

      save_history <- private$save_history

      restart_state <- self$restart_state
      save_restart_step <- private$save_restart_step
      save_restart <- !is.null(restart_state)

      if (save_restart) {
        stop("fixme")
        ## steps_split <- deterministic_steps_restart(steps, save_restart,
        ##                                            private$data)
        ## restart_state <- vector("list", length(save_restart))
        ## y <- vector("list", length(steps_split))
        ## for (i in seq_along(steps_split)) {
        ##   y[[i]] <- model$simulate(steps_split[[i]])
        ##   if (i <= length(save_restart)) {
        ##     s <- model$state()
        ##     restart_state[[i]] <- array_reshape(s, 1, c(nrow(s), 1))
        ##   }
        ## }
        ## y <- array_bind(arrays = y)
        ## restart_state <- array_bind(arrays = restart_state)
      } else {
        if (curr == 0) {
          s <- c(steps[[1]], steps[seq_len(step_index), 2])
        } else {
          s <- steps[curr:step_index, 2]
        }
        y <- model$simulate(s)
        restart_state <- NULL
      }

      if (is.null(index)) {
        y_compare <- y
      } else {
        y_compare <- y[index$run, , , drop = FALSE]
        rownames(y_compare) <- names(index$run)
      }
      log_likelihood <- vnapply(
        seq_along(private$pars), deterministic_likelihood,
        y_compare, private$compare, private$pars, private$data_split)

      if (save_restart) {
        stop("fixme")
      }

      if (save_history) {
        if (is.null(index)) {
          y_history <- y
        } else {
          y_history <- y[index$state, , , drop = FALSE]
          rownames(y_history) <- names(index$state)
        }
        self$history <- list(value = y_history, index = index$predict)
      } else {
        self$history <- NULL
      }

      self$log_likelihood <- log_likelihood
      self$current_step_index <- step_index

      log_likelihood
    },

    fork_multistage = function(pars, transform_state) {
      stop("Writeme")
    }
  ))
