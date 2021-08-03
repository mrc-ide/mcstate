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

    run = function(pars) {
      self$run_many(list(pars))
    },

    run_many = function(pars) {
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
        ## This is a little nastier than expected because we were
        ## never expecting this sort of topology here
        initial_data <- nofilter_initial(pars, private$initial, model$info())
        if (is.list(initial_data)) {
          steps <- particle_steps(steps, initial_data$step)
          model$set_state(initial_data$state, initial_data$step)
        } else {
          model$set_state(initial_data)
        }
      }

      if (is.null(private$index)) {
        index_data <- NULL
      } else {
        index <- nofilter_index(private$index(model$info()))
        model$set_index(index$index)
      }

      y <- model$simulate(steps[, 2], deterministic = TRUE)

      y_compare <- y[index$run, , , drop = FALSE]
      rownames(y_compare) <- rownames(index$data$run)

      ll <- vnapply(seq_len(n_particles), nofilter_likelihood,
                    y_compare, private$compare, pars, private$data_split)

      y_history <- y[index$state, , , drop = FALSE]
      rownames(y_history) <- rownames(index$data$state)

      private$last_model <- model
      private$last_history <- y_history

      ll
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
       state = match(index$state, index_all))
}


nofilter_likelihood <- function(idx, y, compare, pars, data) {
  n_steps <- length(data)
  ll <- numeric(n_steps)
  for (i in seq_len(n_steps)) {
    y_i <- array_drop(y[, idx, i, drop = FALSE], 3L)
    ll[i] <- compare(y_i, data[[i]], pars)
  }
  sum(ll)
}


nofilter_initial <- function(pars, initial, info) {
  init <- lapply(pars, function(p) initial(info(), 1L, p))
  if (length(pars) == 1L) {
    return(init[[1L]])
  }

  ret <- list()
  if ("state" %in% names(init[[1]])) {
    ret$state <- vapply(init, function(x) x$state, init[[1]]$state)
  }

  ret <- list()
  if ("step" %in% names(init[[1]])) {
    ret$step <- vnapply(init, function(x) x$step)
  }

  ret
}
