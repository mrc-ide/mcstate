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

      if (save_history) {
        len <- nrow(steps) + 1L
        state <- model$state(index$predict)
        history_value <- array(NA_real_, c(dim(state), len))
        history_value[, , 1] <- state
        rownames(history_value) <- names(index$predict)
        self$history <- list(
          value = history_value,
          index = index$predict)
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
      private$compare <- compare
      private$save_history <- save_history
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
      self$step(nrow(private$steps))
    },

    step = function(step_index) {
      curr <- self$current_step_index
      check_step(curr, step_index, private$steps)

      model <- self$model
      index <- private$index

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
      log_likelihood <- self$log_likelihood + vnapply(
        seq_along(private$pars), deterministic_likelihood,
        y_compare, private$compare, private$pars,
        private$data_split[idx])

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

    fork_multistage = function(pars, transform_state) {
      stop("Writeme")
    }
  ))


deterministic_likelihood <- function(idx, y, compare, pars, data) {
  ll <- numeric(length(data))
  for (i in seq_along(ll)) {
    y_i <- array_drop(y[, idx, i, drop = FALSE], 3L)
    ll[i] <- compare(y_i, data[[i]], pars[[idx]])
  }
  sum(ll)
}


deterministic_steps_restart <- function(save_restart_step, step_end) {
  i <- match(save_restart_step, step_end)
  j <- which(!is.na(i))

  ## No restart in this block, do the easy exit:
  if (length(j) == 0L) {
    return(list(step_end = step_end, restart = NA_integer_))
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
