##' @title Deterministic particle likelihood
##'
##' @description Create a deterministic version of the
##'   [`mcstate::particle_filter`] object, which runs a single
##'   particle deterministically.
##'
##' @export
particle_deterministic <- R6::R6Class(
  "particle_deterministic",
  cloneable = FALSE,

  public = list(
    ##' @field model The dust model generator being simulated (cannot be
    ##' re-bound)
    model = NULL,

    ##' @description Create the particle filter
    ##'
    ##' @param data The data set to be used for the particle filter,
    ##' created by [particle_filter_data()]. This is essentially
    ##' a [data.frame()] with at least columns `step_start`
    ##' and `step_end`, along with any additional data used in the
    ##' `compare` function, and additional information about how your
    ##' steps relate to time.
    ##'
    ##' @param model A stochastic model to use.  Must be a
    ##' `dust_generator` object.
    ##'
    ##' @param compare A comparison function.  Must take arguments
    ##' `state`, `observed` and `pars` as arguments (though the arguments
    ##' may have different names). `state` is the simulated model state
    ##' (a matrix with as many rows as there are state variables and as
    ##' many columns as there are particles, `data`
    ##' is a `list` of observed data corresponding to the current
    ##' time's row in the `data` object provided here in the
    ##' constructor.  `pars` is any additional parameters passed
    ##' through to the comparison function (via the `pars`
    ##' argument to `$run`). Alternatively, `compare` can be `NULL`
    ##' if your model provides a built-in compile compare function
    ##' (if `model$public_methods$has_compare()` is `TRUE`), which may
    ##' be faster.
    ##'
    ##' @param index An index function. This is used to compute the
    ##' "interesting" indexes of your model. It must be a function of
    ##' one argument, which will be the result of calling the
    ##' `$info()` method on your model. It should return a list
    ##' with elements `run` (indices to return at the end of each
    ##' run, passed through to your compare function) and `state`
    ##' (indices to return if saving state). These indices can overlap
    ##' but do not have to. This argument is optional but using it will
    ##' likely speed up your simulation if you have more than a few
    ##' states as it will reduce the amount of memory copied back and
    ##' forth.
    ##'
    ##' @param initial A function to generate initial conditions. If
    ##' given, then this function must accept 3 arguments: `info`
    ##' (the result of calling `$info()` as for `index`),
    ##' `n_particles` (the number of particles that the particle
    ##' filter is using) and `pars` (parameters passed in in the
    ##' `$run` method via the `pars` argument).  It
    ##' must return a list, which can have the elements `state`
    ##' (initial model state, passed to the particle filter - either a
    ##' vector or a matrix, and overriding the initial conditions
    ##' provided by your model) and `step` (the initial step,
    ##' overriding the first step of your data - this must occur within
    ##' your first epoch in your `data` provided to the
    ##' constructor, i.e., not less than the first element of
    ##' `step_start` and not more than `step_end`). Your function
    ##' can also return a vector or matrix of `state` and not alter
    ##' the starting step, which is equivalent to returning
    ##' `list(state = state, step = NULL)`.
    ##'
    ##' @param n_threads Number of threads to use when running the
    ##' simulation. Defaults to 1, and should not be set higher than the
    ##' number of cores available to the machine. This currently has no
    ##' effect as the simulation will be run in serial on a single
    ##' particle for now.
    initialize = function(data, model, compare,
                          index = NULL, initial = NULL,
                          n_threads = 1L) {
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
      private$compare <- assert_function(compare)
      if (!is.null(index)) {
        private$index <- assert_function(index)
      }
      if (!is.null(initial)) {
        private$initial <- assert_function(initial)
      }
      private$n_threads <- n_threads

      self$model <- model
      lockBinding("model", self)
    },

    ##' @description Run the deterministic particle filter
    ##'
    ##' @param pars A list representing parameters. This will be passed as
    ##' the `pars` argument to your model, to your `compare`
    ##' function, and (if using) to your `initial` function. It must
    ##' be an R list (not vector or `NULL`) because that is what a
    ##' dust model currently requires on initialisation or `$reset` - we
    ##' may relax this later. You may want to put your observation and
    ##' initial parameters under their own keys (e.g.,
    ##' `pars$initial$whatever`), but this is up to you. Extra keys
    ##' are silently ignored by dust models.
    ##'
    ##' @param save_history Logical, indicating if the history of all
    ##' particles should be saved. If saving history, then it can be
    ##' queried later with the `$history` method on the object.
    ##'
    ##' @param save_restart An integer vector of time points to save
    ##' restart infomation for. Not currently supported.
    ##'
    ##' @param min_log_likelihood Not currently supported, exists to match
    ##'   the inteface with [mcstate::particle_filter]. Providing a value
    ##'   larger than -Inf will cause an error.
    ##'
    ##' @return A single numeric value representing the log-likelihood
    ##' (`-Inf` if the model is impossible)
    run = function(pars = list(), save_history = FALSE, save_restart = NULL,
                   min_log_likelihood = -Inf) {
      self$run_many(list(pars), save_history, save_restart, min_log_likelihood)
    },

    ##' @description Run the deterministic particle filter on several
    ##' parameter sets simultaneously. This acts as a wrapper around
    ##' `$run()`, though runs will be in parallel if the object was created
    ##' with `n_threads` greater than 1.
    ##'
    ##' @param pars A list of parameter sets, each element of which would
    ##'  be suitable to pass to `$run()`
    ##'
    ##' @param save_history Logical, indicating if the history of all
    ##' particles should be saved.
    ##'
    ##' @param save_restart An integer vector of time points to save
    ##' restart infomation for. Not currently supported.
    ##'
    ##' @param min_log_likelihood Not currently supported, exists to match
    ##'   the inteface with [mcstate::particle_filter]. Providing a value
    ##'   larger than -Inf will cause an error.
    ##'
    ##' @return A numeric vector of values representing the log-likelihood
    ##' (`-Inf` if the model is impossible), one per parameter set
    run_many = function(pars, save_history = FALSE, save_restart = NULL,
                        min_log_likelihood = -Inf) {
      if (min_log_likelihood > -Inf) {
        stop("'min_log_likelihood' cannot be used with particle_deterministic")
      }
      n_particles <- length(pars)
      steps <- unname(as.matrix(private$data[c("step_start", "step_end")]))
      model <- private$last_model
      resize <- !is.null(model) && n_particles != model$n_particles()

      if (is.null(model) || resize) {
        model <- private$generator$new(pars = pars, step = steps[[1]],
                                       n_particles = NULL,
                                       n_threads = private$n_threads,
                                       seed = NULL,
                                       deterministic = TRUE,
                                       pars_multi = TRUE)
      } else {
        model$update_state(pars = pars, step = steps[[1]])
      }

      if (!is.null(private$initial)) {
        initial_data <- deterministic_initial(pars, private$initial,
                                              model$info())
        steps <- particle_steps(steps, initial_data$step)
        model$update_state(state = initial_data$state,
                           step = initial_data$step)
      }

      if (is.null(private$index)) {
        index <- NULL
      } else {
        ## NOTE: this assumes that all parameterisation result in the
        ## same shape, which is assumed generally.
        index <- deterministic_index(private$index(model$info()[[1L]]))
        model$set_index(index$index)
      }

      if (is.null(save_restart)) {
        y <- model$simulate(c(steps[[1]], steps[, 2]))
        restart_state <- NULL
      } else {
        steps_split <- deterministic_steps_restart(steps, save_restart,
                                                   private$data)
        restart_state <- vector("list", length(save_restart))
        y <- vector("list", length(steps_split))
        for (i in seq_along(steps_split)) {
          y[[i]] <- model$simulate(steps_split[[i]])
          if (i <= length(save_restart)) {
            s <- model$state()
            restart_state[[i]] <- array_reshape(s, 1, c(nrow(s), 1))
          }
        }
        y <- array_bind(arrays = y)
        restart_state <- array_bind(arrays = restart_state)
      }

      if (is.null(index)) {
        y_compare <- y
      } else {
        y_compare <- y[index$run, , , drop = FALSE]
        rownames(y_compare) <- names(index$run)
      }

      ll <- vnapply(seq_len(n_particles), deterministic_likelihood,
                    y_compare, private$compare, pars, private$data_split)

      if (save_history) {
        if (is.null(index)) {
          y_history <- y
        } else {
          y_history <- y[index$state, , , drop = FALSE]
          rownames(y_history) <- names(index$state)
        }
        history <- list(value = y_history, index = index$predict)
      } else {
        history <- NULL
      }

      private$last_model <- model
      private$last_history <- history
      private$last_restart_state <- restart_state

      ll
    },


    ##' @description Extract the current model state, optionally filtering.
    ##' If the model has not yet been run, then this method will throw an
    ##' error. Returns a matrix with the number of rows being the number of
    ##' model states, and the number of columns being the number of
    ##' particles.
    ##'
    ##' @param index_state Optional vector of states to extract
    state = function(index_state = NULL) {
      if (is.null(private$last_model)) {
        stop("Model has not yet been run")
      }
      private$last_model$state(index_state)
    },

    ##' @description Extract the particle trajectories. Requires that
    ##' the model was run with `save_history = TRUE`, which does
    ##' incur a performance cost. This method will throw an error if
    ##' the model has not run, or was run without `save_history =
    ##' TRUE`. Returns a 3d array with dimensions corresponding to (1)
    ##' model state, filtered by `index$run` if provided, (2)
    ##' particle (following `index_particle` if provided), (3)
    ##' time point.
    ##'
    ##' @param index_particle Optional vector of particle indices to return.
    ##' If `NULL` we return all particles' histories.
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

    ##' @description
    ##' Return the full particle filter state at points back in time
    ##' that were saved with the `save_restart` argument to
    ##' `$run()`. If available, this will return a 3d array, with
    ##' dimensions representing (1) particle state, (2) particle index,
    ##' (3) time point. If nested parameters are used then returns a 4d array,
    ##' with dimensions representing (1) particle state, (2) particle index,
    ##' (3) population, (4) time point. This could be quite large, especially
    ##' if you are using the `index` argument to create the particle filter
    ##' and return a subset of all state generally. In the stochastic version,
    ##' this is different the saved trajectories returned by `$history()`
    ##' because earlier saved state is not filtered by later filtering,
    ##' but in the deterministic model we run with a single particle so it
    ##' *is* the same.
    ##'
    ##' @param index_particle Optional vector of particle indices to return.
    ##' If `NULL` we return all particles' states. Practically because the
    ##' only valid value of index_particle is "1", this has no effect and
    ##' it is included primarily for compatibility with the stochastic
    ##' interface.
    restart_state = function(index_particle = NULL) {
      if (is.null(private$last_model)) {
        ## uncovered
        stop("Model has not yet been run")
      }
      restart_state <- private$last_restart_state
      if (is.null(restart_state)) {
        stop("Can't get history as model was run with save_restart = NULL")
      }
      if (!is.null(index_particle)) {
        ## uncovered
        ## NOTE: anything other than 1 here will go poorly; we might
        ## replace with rep(1, length(restart_state)) or at least
        ## validate?
        ##
        ## TODO: nested deterministic filter is not supported so this
        ## is always TRUE; see stochastic filter for the logic when
        ## this is enabled.
        stopifnot(length(dim(restart_state)) == 3)
        restart_state <- restart_state[, index_particle, , drop = FALSE]
      }
      restart_state
    },

    ##' @description
    ##' Return a list of inputs used to configure the deterministic particle
    ##' filter. These correspond directly to the argument names for the
    ##' constructor and are the same as the input arguments.
    inputs = function() {
      list(data = private$data,
           model = private$generator,
           index = private$index,
           initial = private$initial,
           compare = private$compare,
           n_threads = private$n_threads,
           seed = filter_current_seed(private$last_model, NULL))
    },

    ##' @description
    ##' Set the number of threads used by the particle filter (and dust
    ##'   model) after creation. This can be used to allocate additional
    ##'   (or subtract excess) computing power from the deterministic filter
    ##'   Returns (invisibly) the previous value.
    ##'
    ##' @param n_threads The new number of threads to use. You may want to
    ##'   wrap this argument in [dust::dust_openmp_threads()] in order to
    ##'   verify that you can actually use the number of threads
    ##'   requested (based on environment variables and OpenMP support).
    set_n_threads = function(n_threads) {
      prev <- private$n_threads
      private$n_threads <- n_threads
      if (!is.null(private$last_model)) {
        private$last_model$set_n_threads(n_threads)
      }
      invisible(prev)
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
    last_model = NULL,
    last_history = NULL,
    last_restart_state = NULL
  ))


deterministic_index <- function(index) {
  index_all <- union(index$run, index$state)
  list(index = index_all,
       run = set_names(match(index$run, index_all), names(index$run)),
       state = set_names(match(index$state, index_all), names(index$state)),
       predict = index$state)
}


deterministic_likelihood <- function(idx, y, compare, pars, data) {
  n_steps <- length(data)
  ll <- numeric(n_steps)
  for (i in seq_len(n_steps)) {
    y_i <- array_drop(y[, idx, i + 1L, drop = FALSE], 3L)
    ll[i] <- compare(y_i, data[[i]], pars[[idx]])
  }
  sum(ll)
}


deterministic_initial <- function(pars, initial, info) {
  init <- Map(function(p, i) initial(i, 1L, p), pars, info)

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


deterministic_steps_restart <- function(steps, save_restart, data) {
  save_restart_step <- check_save_restart(save_restart, data)
  i <- match(save_restart_step, steps[, 2])
  if (last(i) < nrow(steps)) {
    i <- c(i, nrow(steps))
  }
  j <- rep(seq_along(i), diff(c(0, i)))
  steps_split <- unname(split(steps[, 2], j))
  steps_split[[1]] <- c(steps[[1]], steps_split[[1]])
  steps_split
}
