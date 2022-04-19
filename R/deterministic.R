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

  private = list(
    data = NULL,
    data_split = NULL,
    steps = NULL,
    n_steps = NULL,
    n_threads = NULL,
    initial = NULL,
    index = NULL,
    compare = NULL,
    constant_log_likelihood = NULL,
    last_stages = NULL,
    last_model = NULL,
    last_history = NULL,
    last_state = NULL,
    last_restart_state = NULL
  ),

  public = list(
    ##' @field model The dust model generator being simulated (cannot be
    ##' re-bound)
    model = NULL,

    ##' @field has_multiple_parameters Logical, indicating if the
    ##'   deterministic particle requires multiple parameter sets in a list
    ##'   as inputs, and if it it will produce a vector of likelihoods
    ##'   the same length (read only).  The parameter sets may or may
    ##'   not use the same data (see `has_multiple_data`).
    has_multiple_parameters = NULL,

    ##' @field has_multiple_data Logical, indicating if the deterministic
    ##'   particle simultaneously calculates the likelihood for multiple
    ##'   parameter sets (read only). If `TRUE`, `has_multiple_parameters`
    ##'   will always be `TRUE`.
    has_multiple_data = NULL,

    ##' @field n_parameters The number of parameter sets used by this
    ##'   deterministic particle (read only).  The returned vector of
    ##'   likelihoods will be this length, and if `has_multiple_parameters`
    ##'   is `FALSE` this will be 1.
    n_parameters = NULL,

    ##' @field n_data The number of data sets used by this deterministic
    ##'   particle (read only).  This will either be 1 or the same value as
    ##'   `n_parameters`.
    n_data = NULL,

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
    ##' argument to `$run`).
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
    ##' the process.
    ##'
    ##' @param n_threads Number of threads to use when running the
    ##' simulation. Defaults to 1, and should not be set higher than the
    ##' number of cores available to the machine. This currently has no
    ##' effect as the simulation will be run in serial on a single
    ##' particle for now.
    ##'
    ##' @param n_parameters Number of parameter sets required.  This, along
    ##'   with `data`, controls the interpretation of how the deterministic
    ##'   particle, and importantly will add an additional dimension to
    ##'   most outputs (scalars become vectors, vectors become matrices etc).
    initialize = function(data, model, compare,
                          index = NULL, initial = NULL,
                          constant_log_likelihood = NULL, n_threads = 1L,
                          n_parameters = NULL) {
      if (!is_dust_generator(model)) {
        stop("'model' must be a dust_generator")
      }
      assert_function_or_null(index)
      assert_function_or_null(initial)
      assert_function_or_null(constant_log_likelihood)
      assert_is(data, "particle_filter_data")

      check_compare(compare, model)

      ## NOTE: unlike the particle filter, there is no support for GPU
      ## here (probably never will be)

      self$model <- model
      private$data <- data

      copy_list_and_lock(check_n_parameters(n_parameters, data),
                         self)

      private$steps <- attr(data, "steps")
      private$data_split <- particle_filter_data_split(data, is.null(compare))

      private$compare <- compare
      private$index <- index
      private$initial <- initial
      private$constant_log_likelihood <- constant_log_likelihood

      private$n_threads <- n_threads

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
      filter_run(self, private, pars, save_history, save_restart,
                 min_log_likelihood)
    },

    ##' @description Begin a deterministic run. This is part of the
    ##' "advanced" interface; typically you will want to use `$run()`
    ##' which provides a user-facing wrapper around
    ##' this function. Once created with `$run_begin()`, you should take
    ##' as many steps as needed with `$step()`.
    ##'
    ##' @param pars A list representing parameters. See `$run_many()`
    ##'   for details (and *not* `$run()`)
    ##'
    ##' @param save_history Logical, indicating if the history of all
    ##'   particles should be saved. See `$run()` for details.
    ##'
    ##' @param save_restart Times to save restart state at. See `$run()` for
    ##'   details.
    ##'
    ##' @param min_log_likelihood Not currently supported, exists to match
    ##'   the inteface with [mcstate::particle_filter]. Providing a value
    ##'   larger than -Inf will cause an error.
    ##'
    ##' @return An object of class `particle_deterministic_state`, with methods
    ##' `step` and `end`. This interface is still subject to change.
    run_begin = function(pars, save_history = FALSE, save_restart = NULL,
                         min_log_likelihood = -Inf) {
      assert_scalar_logical(save_history)
      if (min_log_likelihood > -Inf) {
        stop("'min_log_likelihood' cannot be used with particle_deterministic")
      }
      particle_deterministic_state$new(
        pars, self$model, private$last_model[[1]], private$data,
        private$data_split, private$steps, self$has_multiple_parameters,
        private$n_threads, private$initial, private$index, private$compare,
        private$constant_log_likelihood,
        save_history, save_restart)
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
      last(private$last_model)$state(index_state)
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
      state <- private$last_history$value
      if (is.null(state)) {
        stop("Can't get history as model was run with save_history = FALSE")
      }
      if (!is.null(index_particle)) {
        if (length(index_particle) != 1 || index_particle != 1) {
          stop("Invalid value for 'index_particle' may only be 1 (or NULL)")
        }
      }
      state
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
        stop("Model has not yet been run")
      }
      state <- private$last_restart_state
      if (is.null(state)) {
        stop("Can't get history as model was run with save_restart = NULL")
      }
      if (!is.null(index_particle)) {
        if (length(index_particle) != 1 || index_particle != 1) {
          stop("Invalid value for 'index_particle' may only be 1 (or NULL)")
        }
      }
      state
    },

    ##' @description
    ##' Return a list of inputs used to configure the deterministic particle
    ##' filter. These correspond directly to the argument names for the
    ##' constructor and are the same as the input arguments.
    inputs = function() {
      list(data = private$data,
           model = self$model,
           index = private$index,
           initial = private$initial,
           compare = private$compare,
           constant_log_likelihood = private$constant_log_likelihood,
           n_threads = private$n_threads,
           seed = filter_current_seed(last(private$last_model), NULL))
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
      particle_filter_set_n_threads(private, n_threads)
    }
  ))
