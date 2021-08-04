particle_nofilter <- R6::R6Class(
  "particle_nofilter",
  cloneable = FALSE,

  public = list(
    model = NULL,
    initialize = function(data, model, compare,
                          index = NULL, initial = NULL,
                          n_threads = 1L, seed = NULL) {
      if (!is_dust_generator(model)) {
        stop("'model' must be a dust_generator")
      }
      if (!is.null(index) && !is.function(index)) {
        stop("'index' must be function if not NULL")
      }
      if (!is.null(initial) && !is.function(initial)) {
        stop("'initial' must be function if not NULL")
      }
      assert_is(data, "particle_filter_data")

      if (inherits(data, "particle_filter_data_nested")) {
        stop("nested mode not yet supported")
      }

      private$generator <- model
      private$data <- data
      private$data_split <- df_to_list_of_lists(data)
      private$compare <- assert_is(compare, "function")
      if (!is.null(index)) {
        private$index <- assert_is(index, "function")
      }
      if (!is.null(initial)) {
        private$initial <- assert_is(initial, "function")
      }
      private$n_threads <- n_threads
      private$seed <- seed
    },

    run = function(pars = list(), save_history = FALSE, save_restart = NULL) {
      self$run_many(list(pars), save_history, save_restart)
    },

    run_many = function(pars, save_history = FALSE, save_restart = NULL) {
      if (!is.null(save_restart)) {
        stop("'save_restart' cannot be used with the deterministic nofilter")
      }
      n_particles <- length(pars)
      steps <- unname(as.matrix(private$data[c("step_start", "step_end")]))
      model <- private$last_model
      resize <- !is.null(model) && n_particles != model$n_particles()

      if (is.null(model) || resize) {
        seed <- if (resize) model$rng_state() else private$seed
        model <- private$generator$new(pars = pars, step = steps[[1]],
                                       n_particles = NULL,
                                       n_threads = private$n_threads,
                                       seed = private$seed,
                                       pars_multi = TRUE)
      } else {
        model$reset(pars, steps[[1]])
      }

      if (!is.null(private$initial)) {
        initial_data <- nofilter_initial(pars, private$initial, model$info())
        steps <- particle_steps(steps, initial_data$step)
        model$set_state(initial_data$state, initial_data$step,
                        deterministic = TRUE)
      }

      if (is.null(private$index)) {
        index <- NULL
      } else {
        index <- nofilter_index(private$index(model$info()))
        model$set_index(index$index)
      }

      y <- model$simulate(c(steps[[1]], steps[, 2]), deterministic = TRUE)

      if (is.null(index)) {
        y_compare <- y
      } else {
        y_compare <- y[index$run, , , drop = FALSE]
        rownames(y_compare) <- rownames(index$data$run)
      }

      ll <- vnapply(seq_len(n_particles), nofilter_likelihood,
                    y_compare, private$compare, pars, private$data_split)

      if (save_history) {
        if (is.null(index)) {
          y_history <- y
        } else {
          y_history <- y[index$state, , , drop = FALSE]
          rownames(y_history) <- rownames(index$data$state)
        }
        history <- list(value = y_history, index = index$predict)
      } else {
        history <- NULL
      }

      private$last_model <- model
      private$last_history <- history

      ll
    },

    state = function(index_state = NULL) {
      if (is.null(private$last_model)) {
        stop("Model has not yet been run")
      }
      private$last_model$state(index_state)
    },

    history = function(index_particle = NULL) {
      if (is.null(private$last_model)) {
        stop("Model has not yet been run")
      }
      if (is.null(private$last_history)) {
        stop("Can't get history as model was run with save_history = FALSE")
      }
      if (is.null(index_particle)) {
        private$last_history$value
      } else {
        array_drop(
          private$last_history$value[, index_particle, , drop = FALSE], 2L)
      }
    },

    inputs = function() {
      if (is.null(private$last_model)) {
        seed <- private$seed
      } else {
        seed <- private$last_model$rng_state(first_only = TRUE)
      }
      list(data = private$data,
           model = self$model,
           n_particles = self$n_particles,
           index = private$index,
           initial = private$initial,
           compare = private$compare,
           device_config = private$device_config,
           n_threads = private$n_threads,
           seed = seed)
    }
  ),
  private = list(
    generator = NULL,
    data = NULL,
    data_split = NULL,
    steps = NULL,
    n_steps = NULL,
    n_threads = NULL,
    initial = NULL,
    index = NULL,
    compare = NULL,
    seed = NULL,
    last_model = NULL,
    last_history = NULL
  ))


nofilter_index <- function(index) {
  index_all <- union(index$run, index$state)
  list(index = index_all,
       run = match(index$run, index_all),
       state = match(index$state, index_all),
       predict = index$state)
}


nofilter_likelihood <- function(idx, y, compare, pars, data) {
  n_steps <- length(data)
  ll <- numeric(n_steps)
  for (i in seq_len(n_steps)) {
    y_i <- array_drop(y[, idx, i + 1L, drop = FALSE], 3L)
    ll[i] <- compare(y_i, data[[i]], pars)
  }
  sum(ll)
}


nofilter_initial <- function(pars, initial, info) {
  init <- lapply(pars, function(p) initial(info(), 1L, p))

  ret <- list()
  if (is.list(init[[1]])) {
    if ("state" %in% names(init[[1]])) {
      ret$state <- vapply(init, function(x) x$state, init[[1]]$state)
    }

    if ("step" %in% names(init[[1]])) {
      ret$step <- vnapply(init, function(x) x$step)
    }
  } else {
    ret$state <- vapply(init, identity, init[[1]])
  }

  ret
}
