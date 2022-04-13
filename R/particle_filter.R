##' @title Particle filter
##'
##' @description Create a `particle_filter` object for running
##'   and interacting with a particle filter.  A higher-level
##'   interface will be implemented later.
##'
##' @export
##' @importFrom R6 R6Class
##' @examples
##' # A basic SIR model included in the dust package
##' gen <- dust::dust_example("sir")
##'
##' # Some data that we will fit to, using 1 particle:
##' sir <- gen$new(pars = list(), step = 0, n_particles = 1)
##' dt <- 1 / 4
##' day <- seq(1, 100)
##' incidence <- rep(NA, length(day))
##' true_history <- array(NA_real_, c(5, 1, 101))
##' true_history[, 1, 1] <- sir$state()
##' for (i in day) {
##'   state_start <- sir$state()
##'   sir$run(i / dt)
##'   state_end <- sir$state()
##'   true_history[, 1, i + 1] <- state_end
##'   # Reduction in S
##'   incidence[i] <- state_start[1, 1] - state_end[1, 1]
##' }
##'
##' # Convert this into our required format:
##' data_raw <- data.frame(day = day, incidence = incidence)
##' data <- particle_filter_data(data_raw, "day", 4)
##'
##' # A comparison function
##' compare <- function(state, observed, pars = NULL) {
##'   if (is.null(pars$exp_noise)) {
##'     exp_noise <- 1e6
##'   } else {
##'     exp_noise <- pars$exp_noise
##'   }
##'   incidence_modelled <- state[1,]
##'   incidence_observed <- observed$incidence
##'   lambda <- incidence_modelled +
##'     rexp(length(incidence_modelled), exp_noise)
##'   dpois(incidence_observed, lambda, log = TRUE)
##' }
##'
##' # Construct the particle_filter object with 100 particles
##' p <- particle_filter$new(data, gen, 100, compare)
##' p$run(save_history = TRUE)
##'
##' # Our simulated trajectories, with the "real" data superimposed
##' history <- p$history()
##' matplot(data_raw$day, t(history[1, , -1]), type = "l",
##'         xlab = "Time", ylab = "State",
##'         col = "#ff000022", lty = 1, ylim = range(history))
##' matlines(data_raw$day, t(history[2, , -1]), col = "#ffff0022", lty = 1)
##' matlines(data_raw$day, t(history[3, , -1]), col = "#0000ff22", lty = 1)
##' matpoints(data_raw$day, t(true_history[1:3, , -1]), pch = 19,
##'           col = c("red", "yellow", "blue"))
particle_filter <- R6::R6Class(
  "particle_filter",
  cloneable = FALSE,

  private = list(
    ## Control over the data
    data = NULL,
    data_split = NULL,
    steps = NULL,
    ## Functions used for initial conditions, data comparisons and indices
    index = NULL,
    initial = NULL,
    compare = NULL,
    constant_log_likelihood = NULL,
    gpu_config = NULL,
    ## Control for dust
    seed = NULL,
    n_threads = NULL,
    ## Updated when the model is run
    last_stages = NULL,
    last_model = NULL,
    last_state = NULL,
    last_history = NULL,
    last_restart_state = NULL
  ),

  public = list(
    ##' @field model The dust model generator being simulated (cannot be
    ##' re-bound)
    model = NULL,

    ##' @field n_particles Number of particles used (read only)
    n_particles = NULL,

    ## We don't actually need to expose most of these I think...
    has_multiple_data = NULL,
    has_multiple_parameters = NULL,
    n_data = NULL,
    n_parameters = NULL,

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
    ##' @param n_particles The number of particles to simulate
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
    ##' @param constant_log_likelihood An optional function, taking the
    ##' model parameters, that computes the constant part of the
    ##' log-likelihood value (if any).  You can use this where your
    ##' likelihood depends both on the time series (via `data`) but also
    ##' on some non-temporal data.  You should bind any non-parameter
    ##' dependencies into this closure.  This is applied at the
    ##' beginning of the filter run, so represents the initial
    ##' condition of the marginal log likelihood value propagated by
    ##' the filter.
    ##'
    ##' @param n_threads Number of threads to use when running the
    ##' simulation. Defaults to 1, and should not be set higher than the
    ##' number of cores available to the machine.
    ##'
    ##' @param n_parameters Number of parameter sets required.  This, along
    ##'   with `data`, controls the interpretation of how the particle
    ##'   filter, and importantly will add an additional dimension to
    ##'   most outputs (scalars become vectors, vectors become matrices etc).
    ##'
    ##' @param seed Seed for the random number generator on initial
    ##' creation. Can be `NULL` (to initialise using R's random number
    ##' generator), a positive integer, or a raw vector - see [`dust::dust`]
    ##' and [`dust::dust_rng`] for more details. Note that the random number
    ##' stream is unrelated from R's random number generator, except for
    ##' initialisation with `seed = NULL`.
    ##'
    ##' @param gpu_config GPU configuration, typically an integer
    ##' indicating the device to use, where the model has GPU support.
    ##' An error is thrown if the device id given is larger than those
    ##' reported to be available (note that CUDA numbers devices from 0,
    ##' so that '0' is the first device, so on). See the method `$gpu_info()`
    ##' for available device ids; this can be called before object creation
    ##' as `model$public_methods$gpu_info()`.
    ##' For additional control, provide a list with elements `device_id`
    ##' and `run_block_size`. Further options (and validation) of this
    ##' list will be added in a future version!
    initialize = function(data, model, n_particles, compare,
                          index = NULL, initial = NULL,
                          constant_log_likelihood = NULL,
                          n_threads = 1L, seed = NULL,
                          n_parameters = NULL,
                          gpu_config = NULL) {
      if (!is_dust_generator(model)) {
        stop("'model' must be a dust_generator")
      }
      assert_function_or_null(index)
      assert_function_or_null(initial)
      assert_function_or_null(constant_log_likelihood)
      assert_is(data, "particle_filter_data")

      check_compare(compare, model)

      if (!is.null(gpu_config)) {
        if (!model$public_methods$has_gpu_support(TRUE)) {
          stop(paste("'gpu_config' provided, but 'model' does not have",
                     "GPU support"))
        }
      }

      self$model <- model
      private$data <- data

      copy_list_and_lock(check_n_parameters(n_parameters, data),
                         self)

      private$steps <- attr(data, "steps")
      private$data_split <- particle_filter_data_split(data, is.null(compare),
                                                       self$n_parameters)

      private$compare <- compare
      private$gpu_config <- gpu_config
      private$index <- index
      private$initial <- initial
      private$constant_log_likelihood <- constant_log_likelihood

      self$n_particles <- assert_scalar_positive_integer(n_particles)
      private$n_threads <- assert_scalar_positive_integer(n_threads)
      private$seed <- seed

      lockBinding("model", self)
      lockBinding("n_particles", self)
    },

    ##' @description Run the particle filter
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
    ##' restart infomation for. These are in terms of your underlying time
    ##' variable (the `time` column in [particle_filter_data()]) not in
    ##' terms of steps. The state will be saved after the particle
    ##' filtering operation (i.e., at the end of the step).
    ##'
    ##' @param min_log_likelihood Optionally, a numeric value representing the
    ##' smallest likelihood we are interested in. If given and the particle
    ##' filter drops below this number, then we terminate early and return
    ##' `-Inf`. In this case, history and final state cannot be returned
    ##' from the filter. This is primarily intended for use with
    ##' [mcstate::pmcmc] where we can avoid computing likelihoods that
    ##' will certainly be rejected. Only suitable for use where
    ##' log-likelihood increments (with the `compare` function) are always
    ##' negative. This is the case if you use a normalised discrete
    ##' distribution, but not necessarily otherwise. If using a nested
    ##' filter this can be a single number (in which case the exit is
    ##' when the sum of log-likelihoods drops below this threshhold) or
    ##' a vector of numbers the same length as `pars` (in which case exit
    ##' occurs when all numbers drop below this threshhold).
    ##'
    ##' @return A single numeric value representing the log-likelihood
    ##' (`-Inf` if the model is impossible)
    run = function(pars = list(), save_history = FALSE, save_restart = NULL,
                   min_log_likelihood = NULL) {
      filter_run(self, private, pars, save_history, save_restart,
                 min_log_likelihood)
    },

    ##' @description Begin a particle filter run. This is part of the
    ##' "advanced" interface for the particle filter; typically you will
    ##' want to use `$run()` which provides a user-facing wrapper around
    ##' this function. Once created with `$run_begin()`, you should take
    ##' as many steps as needed with `$step()`.
    ##'
    ##' @param pars A list representing parameters. See `$run()` for details.
    ##'
    ##' @param save_history Logical, indicating if the history of all
    ##' particles should be saved. See `$run()` for details.
    ##'
    ##' @param save_restart Times to save restart state at. See `$run()` for
    ##' details.
    ##'
    ##' @param min_log_likelihood Optionally, a numeric value representing the
    ##' smallest likelihood we are interested in. See `$run()` for details.
    ##'
    ##' @return An object of class `particle_filter_state`, with methods
    ##' `step` and `end`. This interface is still subject to change.
    run_begin = function(pars = list(), save_history = FALSE,
                         save_restart = NULL, min_log_likelihood = NULL) {
      assert_scalar_logical(save_history)
      min_log_likelihood <- min_log_likelihood %||% -Inf
      particle_filter_state$new(
        pars, self$model, private$last_model[[1]], private$data,
        private$data_split, private$steps, self$n_particles,
        self$has_multiple_parameters, private$n_threads,
        private$initial, private$index, private$compare,
        private$constant_log_likelihood, private$gpu_config, private$seed,
        min_log_likelihood, save_history, save_restart)
    },

    ##' @description Extract the current model state, optionally filtering.
    ##' If the model has not yet been run, then this method will throw an
    ##' error. Returns a matrix with the number of rows being the number of
    ##' model states, and the number of columns being the number of
    ##' particles.
    ##'
    ##' @param index_state Optional vector of states to extract
    state = function(index_state = NULL) {
      if (is.null(private$last_state)) {
        stop("Model has not yet been run")
      }
      ## TODO (#173): should get an option to take a single trajectory
      private$last_state(index_state)
    },

    ##' @description Extract the particle trajectories. Requires that
    ##' the model was run with `save_history = TRUE`, which does
    ##' incur a performance cost. This method will throw an error if
    ##' the model has not run, or was run without `save_history =
    ##' TRUE`. Returns a 3d array with dimensions corresponding to (1)
    ##' model state, filtered by `index$run` if provided, (2)
    ##' particle (following `index_particle` if provided), (3)
    ##' time point. If nested parameters used then returns a 4d array with
    ##' dimensions corresponding to (1) model state, (2) particle, (3)
    ##' population, (4) time point.
    ##'
    ##' @param index_particle Optional vector of particle indices to return.
    ##' If nested parameters used then a vector will be replicated to a matrix
    ##' with number of columns equal to number of populations, otherwise a
    ##' matrix can be supplied.
    ##' If `NULL` we return all particles' histories.
    history = function(index_particle = NULL) {
      if (is.null(private$last_model)) {
        stop("Model has not yet been run")
      }
      if (is.null(private$last_history)) {
        stop("Can't get history as model was run with save_history = FALSE")
      }

      history_value <- private$last_history$value
      history_order <- private$last_history$order
      history_index <- private$last_history$index

      ny <- nrow(history_value)

      if (length(dim(history_value)) == 4) {
        history_nested(history_value, history_order, history_index,
                       index_particle)
      } else {
        history_single(history_value, history_order, history_index,
                       index_particle)
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
    ##' and return a subset of all state generally. It is also
    ##' different to the saved trajectories returned by `$history()`
    ##' because earlier saved state is not filtered by later filtering
    ##' (in the history we return the tree of history representing the
    ##' histories of the _final_ particles, here we are returning all
    ##' particles at the requested point, regardless if they appear in
    ##' the set of particles that make it to the end of the
    ##' simulation).
    ##'
    ##' @param index_particle Optional vector of particle indices to return.
    ##' If `NULL` we return all particles' states.
    restart_state = function(index_particle = NULL) {
      if (is.null(private$last_model)) {
        stop("Model has not yet been run")
      }
      restart_state <- private$last_restart_state
      if (is.null(restart_state)) {
        stop("Can't get history as model was run with save_restart = NULL")
      }
      if (!is.null(index_particle)) {
        if (length(dim(restart_state)) == 4) {
          restart_state <- restart_state[, index_particle, , , drop = FALSE]
        } else {
          restart_state <- restart_state[, index_particle, , drop = FALSE]
        }
      }
      restart_state
    },

    ##' @description
    ##' Return a list of inputs used to configure the particle
    ##' filter. These correspond directly to the argument names for the
    ##' particle filter constructor and are the same as the input
    ##' argument with the exception of `seed`, which is the state of
    ##' the rng if it has been used (this can be used as a seed to
    ##' restart the model).
    inputs = function() {
      if (self$has_multiple_parameters) {
        n_parameters <- self$n_parameters
      } else {
        n_parameters <- NULL
      }
      list(data = private$data,
           model = self$model,
           n_particles = self$n_particles,
           index = private$index,
           initial = private$initial,
           compare = private$compare,
           constant_log_likelihood = private$constant_log_likelihood,
           gpu_config = private$gpu_config,
           n_threads = private$n_threads,
           n_parameters = n_parameters,
           seed = filter_current_seed(last(private$last_model), private$seed))
    },

    ##' @description
    ##' Set the number of threads used by the particle filter (and dust
    ##'   model) after creation. This can be used to allocate additional
    ##'   (or subtract excess) computing power from a particle filter.
    ##'   Returns (invisibly) the previous value.
    ##'
    ##' @param n_threads The new number of threads to use. You may want to
    ##'   wrap this argument in [dust::dust_openmp_threads()] in order to
    ##'   verify that you can actually use the number of threads
    ##'   requested (based on environment variables and OpenMP support).
    set_n_threads = function(n_threads) {
      particle_filter_set_n_threads(private, n_threads)
    }
  ))


##' @importFrom stats runif
particle_resample <- function(weights) {
  if (is.matrix(weights)) {
    return(apply(weights, 2, particle_resample))
  }
  n <- length(weights)
  u <- runif(1, 0, 1 / n) + seq(0, by = 1 / n, length.out = n)
  cum_weights <- cumsum(weights / sum(weights))
  findInterval(u, cum_weights) + 1L
}


## Private helper for reconstituting a particle filter from its
## `$inputs()` data, but possibly changing the seed
particle_filter_from_inputs <- function(inputs, seed = NULL) {
  if (is.null(inputs$n_particles)) {
    particle_deterministic$new(
      data = inputs$data,
      model = inputs$model,
      compare = inputs$compare,
      index = inputs$index,
      initial = inputs$initial,
      constant_log_likelihood = inputs$constant_log_likelihood,
      n_threads = inputs$n_threads)
  } else {
    particle_filter$new(
      data = inputs$data,
      model = inputs$model,
      n_particles = inputs$n_particles,
      compare = inputs$compare,
      gpu_config = inputs$gpu_config,
      index = inputs$index,
      initial = inputs$initial,
      constant_log_likelihood = inputs$constant_log_likelihood,
      n_threads = inputs$n_threads,
      seed = seed %||% inputs$seed)
  }
}


scale_log_weights <- function(log_weights) {
  log_weights[is.nan(log_weights)] <- -Inf
  max_log_weights <- max(log_weights)
  if (!is.finite(max_log_weights)) {
    ## if all log_weights at a time-step are -Inf, this should
    ## terminate the particle filter and output the marginal
    ## likelihood estimate as -Inf
    average <- -Inf
    weights <- rep(NaN, length(log_weights))
  } else {
    ## calculation of weights, there is some rescaling here to avoid
    ## issues where exp(log_weights) might give computationally zero
    ## values
    weights <- exp(log_weights - max_log_weights)
    average <- log(mean(weights)) + max_log_weights
  }
  list(weights = weights, average = average)
}


particle_steps <- function(steps, step_start) {
  if (!is.null(step_start)) {
    assert_integer(step_start)
    if (min(step_start) < steps[1, 1, drop = TRUE]) {
      stop(sprintf(
        "'step_start' must be >= %d (the first value of data$step_start)",
        steps[1, 1, drop = TRUE]))
    }
    if (max(step_start) > steps[1, 2, drop = TRUE]) {
      stop(sprintf(
        "'step_start' must be <= %d (the first value of data$step_end)",
        steps[1, 2, drop = TRUE]))
    }
    steps[1, 1] <- max(step_start)
  }
  steps
}


is_dust_generator <- function(x) {
  inherits(x, "R6ClassGenerator") &&
    identical(attr(x, which = "name", exact = TRUE), "dust_generator")
}


check_save_restart <- function(save_restart, data) {

  if (is.null(save_restart)) {
    return(integer(0))
  }
  assert_strictly_increasing(save_restart)
  assert_is(data, "particle_filter_data")

  time_end <- attr(data, "times")[, 2]
  i <- match(save_restart, time_end)
  if (anyNA(i)) {
    stop(sprintf("'save_restart' contains times not in '%s': %s",
                 attr(data, "time"),
                 paste(save_restart[is.na(i)], collapse = ", ")))
  }

  data$step_end[i]
}


history_single <- function(history_value, history_order, history_index,
                           index_particle) {
  ny <- nrow(history_value)

  if (is.null(history_order)) {
    if (is.null(index_particle)) {
      ret <- history_value
    } else {
      ret <- history_value[, index_particle, , drop = FALSE]
    }
  } else {
    if (is.null(index_particle)) {
      index_particle <- seq_len(ncol(history_value))
    }

    np <- length(index_particle)
    nt <- ncol(history_order)

    idx <- matrix(NA_integer_, np, nt)
    for (i in rev(seq_len(ncol(idx)))) {
      index_particle <- idx[, i] <- history_order[index_particle, i]
    }

    cidx <- cbind(seq_len(ny),
                  rep(idx, each = ny),
                  rep(seq_len(nt), each = ny * np))
    ret <- array(history_value[cidx], c(ny, np, nt))
  }
  rownames(ret) <- names(history_index)
  ret
}

## This function handles the nested/non-nested case but also the
## compiled/non-compiled case. In the compiled case we already have
## our history nicely ordered (that is the states convered into a tree
## based on the history of particle sampling) and history_order is
## NULL.
history_nested <- function(history_value, history_order, history_index,
                           index_particle) {
  ny <- nrow(history_value)
  npop <- nlayer(history_value)

  if (is.null(history_order)) {
    ## Compiled particle filter; no ordering needed (or available)
    if (is.null(index_particle)) {
      ret <- history_value
    } else if (!is.matrix(index_particle)) {
      ret <- history_value[, index_particle, , , drop = FALSE]
    } else {
      if (!ncol(index_particle) == npop) {
        stop(sprintf("'index_particle' should have %d columns", npop))
      }
      d <- dim(history_value)
      d[[2L]] <- nrow(index_particle)
      ret <- array(NA_real_, d)
      for (i in seq_len(npop)) {
        ret[, , i, ] <- history_value[, index_particle[, i], i, ]
      }
    }
  } else {
    ## mcstate particle filter; need to sort the history
    nt <- nlayer(history_order)

    if (is.null(index_particle)) {
      index_particle <- matrix(seq_len(ncol(history_value)),
                               ncol(history_value), npop)
    } else {
      if (is.matrix(index_particle)) {
        if (!ncol(index_particle) == npop) {
          stop(sprintf("'index_particle' should have %d columns", npop))
        }
      } else {
        index_particle <- matrix(index_particle,
                                 nrow = length(index_particle),
                                 ncol = npop)
      }
    }

    np <- nrow(index_particle)

    idx <- array(NA_integer_, c(np, npop, nt))
    for (i in rev(seq_len(nlayer(idx)))) {
      for (j in seq_len(npop)) {
        idx[, j, i] <- history_order[, j, i][index_particle[, j]]
      }
      index_particle <- matrix(idx[, , i], nrow = np, ncol = npop)
    }

    ret <- array(NA, c(ny, np, npop, nt))
    for (i in seq_len(npop)) {
      cidx <- cbind(seq_len(ny),
                    rep(idx[, i, ], each = ny),
                    rep(seq_len(nt), each = ny * np))
      ret[, , i, ] <- history_value[, , i, ][cidx]
    }
  }
  rownames(ret) <- names(history_index)
  ret
}


filter_current_seed <- function(model, seed) {
  if (!is.null(model)) {
    seed <- model$rng_state(first_only = TRUE)
  }
  seed
}


filter_run <- function(self, private, pars, save_history, save_restart,
                       min_log_likelihood) {
  assert_scalar_logical(save_history)
  if (self$has_multiple_parameters) {
    n_parameters <- self$n_parameters
    pars <- particle_filter_pars_nested(pars, n_parameters) # TODO: rename
  }
  private$last_stages <-
    particle_filter_check_multistage_pars(pars, private$last_stages)
  if (inherits(pars, "multistage_parameters")) {
    filter_run_multistage(self, private, pars, save_history, save_restart,
                          min_log_likelihood)
  } else {
    filter_run_simple(self, private, pars, save_history, save_restart,
                      min_log_likelihood)
  }
}


filter_run_simple <- function(self, private, pars,
                              save_history, save_restart,
                              min_log_likelihood) {
  obj <- self$run_begin(pars, save_history, save_restart,
                        min_log_likelihood = min_log_likelihood)
  obj$run()
  private$last_history <- obj$history
  private$last_model <- list(obj$model)
  private$last_state <- function(index) obj$model$state(index)
  private$last_restart_state <- obj$restart_state
  obj$log_likelihood
}


filter_run_multistage <- function(self, private, pars,
                                  save_history, save_restart,
                                  min_log_likelihood) {
  stages <- filter_check_times(pars, private$data, save_restart)

  models <- private$last_model %||% vector("list", length(stages))
  history <- vector("list", length(stages))
  restart <- vector("list", length(stages))

  for (i in seq_along(stages)) {
    if (i == 1) {
      obj <- self$run_begin(
        stages[[i]]$pars, save_history, save_restart, min_log_likelihood)
    } else {
      obj <- obj$fork_multistage(
        models[[i]], stages[[i]]$pars, stages[[i]]$transform_state)
    }
    obj$step(stages[[i]]$step_index)
    models[[i]] <- obj$model
    history[i] <- list(obj$history)
    restart[i] <- list(obj$restart_state)
  }

  ## Push the final rng state into the first version of the model,
  ## completing the cycle.
  models[[1]]$set_rng_state(last(models)$rng_state())

  ## We return this first model in the sequence as that's where
  ## the next run will start from, but state from the last model
  ## because that's where we got to.
  private$last_model <- models
  private$last_state <- function(index) last(models)$state(index)

  if (save_history) {
    private$last_history <- join_histories(history, stages)
  } else {
    private$last_history <- NULL
  }

  if (!is.null(save_restart)) {
    private$last_restart_state <- join_restart_state(restart, stages)
  } else {
    private$last_restart_state <- NULL
  }

  obj$log_likelihood
}


## There are several bits of cleanup that need to happen for the
## parameters in nested case:
##
## * validate we have an unnamed list of the correct length
## * if multistage, then invert the nesting to convert from a list of
##   multistage parameter objects into a multistage parameter of lists
particle_filter_pars_nested <- function(pars, n_populations) {
  if (!is.null(names(pars))) {
    stop("Expected an unnamed list of parameters for 'pars'")
  }
  if (length(pars) != n_populations) {
    stop(sprintf("'pars' must have length %d", n_populations))
  }

  is_multistage <- vlapply(pars, inherits, "multistage_parameters")
  if (!any(is_multistage)) {
    return(pars)
  }
  if (any(!is_multistage)) {
    stop("'pars' must be either all multistage or all non-multistage")
  }

  ret <- pars[[1L]]
  if (length(pars) > 1 && any(lengths(pars) != length(ret))) {
    stop(sprintf(
      "Incompatible numbers of stages in pars: found %s stages",
      paste(sort(unique(lengths(pars))), collapse = ", ")))
  }

  for (i in seq_along(ret)) {
    if (i > 1 && length(pars) > 1) {
      err_start <- vlapply(pars[-1], function(x)
        x[[i]]$start != ret[[i]]$start)
      if (any(err_start)) {
        stop(sprintf("Incompatible 'start' time at phase %d", i))
      }
      err_transform <- vlapply(pars[-1], function(x)
        !identical(x[[i]]$transform_state, ret[[i]]$transform_state))
      if (any(err_transform)) {
        stop(sprintf("Incompatible 'transform_state' at phase %d", i))
      }
    }

    p <- lapply(pars, function(x) x[[i]]$pars)
    is_null <- vlapply(p, is.null)
    err_pars <- is_null[-1] != is_null[[1]]
    if (any(err_pars)) {
      stop(sprintf("Incompatible 'pars' at phase %d", i))
    }
    if (!all(is_null)) {
      ret[[i]]$pars <- p
    }
  }

  ret
}


particle_filter_set_n_threads <- function(private, n_threads) {
  prev <- private$n_threads
  private$n_threads <- n_threads
  for (m in private$last_model) {
    if (!is.null(m)) {
      m$set_n_threads(n_threads)
    }
  }
  invisible(prev)
}


particle_filter_check_multistage_pars <- function(pars, n_stages_prev) {
  is_multistage <- inherits(pars, "multistage_parameters")
  n_stages_given <- if (is_multistage) length(pars) else 1L
  if (!is.null(n_stages_prev) && n_stages_prev != n_stages_given) {
    if (n_stages_prev == 1) {
      stop(sprintf(
        "Expected single-stage parameters (but given one with %d stages)",
        n_stages_given))
    } else {
      stop(sprintf(
        "Expected multistage_pars with %d stages (but given one with %d)",
        n_stages_prev, n_stages_given))
    }
  }
  n_stages_given
}


check_compare <- function(compare, model) {
  if (is.null(compare)) {
    if (!model$public_methods$has_compare()) {
      stop("Your model does not have a built-in 'compare' function")
    }
  } else {
    assert_function(compare)
  }
}


check_n_parameters <- function(n_parameters, data) {
  has_multiple_data <- inherits(data, "particle_filter_data_nested")
  if (has_multiple_data) {
    n_data <- length(attr(data, "populations"))
  } else {
    n_data <- 1L
  }

  if (is.null(n_parameters)) {
    has_multiple_parameters <- has_multiple_data
    n_parameters <- n_data
  } else {
    assert_scalar_positive_integer(n_parameters)
    if (has_multiple_data && n_parameters != n_data) {
      stop(paste("To match the number of populations in your data,",
                 sprintf("n_parameters must be %d (if not NULL)", n_data)))
    }
    has_multiple_parameters <- TRUE
    n_parameters <- n_parameters
  }
  list(has_multiple_parameters = has_multiple_parameters,
       has_multiple_data = has_multiple_data,
       n_parameters = n_parameters,
       n_data = n_data)
}
