##' Describe a single parameter for use within IF2. Note that
##' the name is not set here, but will end up being naturally defined
##' when used with [`mcstate::if2_parameters`], which collects
##' these together for use with [mcstate::if2()].
##'
##' @title Describe single IF2 parameter
##'
##' @param name Name for the parameter (a string)
##'
##' @param initial Initial value of the parameter
##'
##' @param min Optional minimum value for the parameter (otherwise
##'   `-Inf`). If given, then `initial` must be at least this
##'   value.
##'
##' @param max Optional max value for the parameter (otherwise
##'   `Inf`). If given, then `initial` must be at most this
##'   value.
##'
##' @param discrete Deprecated. Logical, indicating if this parameter is
##'   discrete. If `TRUE` then the parameter will be rounded
##'   after a new parameter is proposed.
##'
##' @param integer Logical, indicating if this parameter is
##'   integer. If `TRUE` then the parameter will be rounded
##'   after a new parameter is proposed.
##'
##' @param prior A prior function (if not given an improper flat prior
##'   is used - be careful!). It must be a function that takes a
##'   single argument, being the value of this parameter. If given,
##'   then `prior(initial)` must evaluate to a finite value.
##'
##' @export
##' @examples
##' mcstate::if2_parameter("a", 0.1)
if2_parameter <- function(name, initial,
                          min = -Inf, max = Inf, discrete,
                          integer = FALSE,
                          prior = NULL) {
  if (!missing(discrete)) {
    .Deprecated("integer", old = "discrete")
    integer <- discrete
  }
  assert_scalar_character(name)
  assert_scalar_logical(integer)

  if (initial < min) {
    stop(sprintf("'initial' must be >= 'min' (%s)", min))
  }
  if (initial > max) {
    stop(sprintf("'initial' must be <= 'max' (%s)", max))
  }
  if (is.null(prior)) {
    prior <- function(p) 0
  } else {
    assert_is(prior, "function")
    value <- tryCatch(
      prior(initial),
      error = function(e)
        stop(sprintf(
          "Prior function for '%s' failed to evaluate initial value: %s",
          name, e$message)))
    if (!is.finite(value)) {
      stop(sprintf(
        "Prior function for '%s' returned a non-finite value on initial value",
        name))
    }
  }

  ret <- list(name = name, initial = initial, min = min, max = max,
              integer = integer, prior = prior)
  class(ret) <- "if2_parameter"
  ret
}


##' @title if2_parameters
##'
##' @description Construct parameters for use with
##'   [mcstate::if2()]. This creates a utility object that is used
##'   internally to work with parameters. Most users only need to
##'   construct this object, but see the examples for how it can be
##'   used.
##'
##' @export
##' @examples
##' # Construct an object with two parameters:
##' pars <- mcstate::if2_parameters$new(
##'   list(mcstate::if2_parameter("a", 0.1, min = 0, max = 1,
##'                                 prior = function(a) log(a)),
##'        mcstate::if2_parameter("b", 0, prior = dnorm)))
##'
##' # Initial parameters
##' pars$initial()
##'
##' # Create the initial parameter set
##' n_par_sets <- 5
##' pars_sd <- list("a" = 0.02, "b" = 0.02)
##' p_mat <- pars$walk_initialise(n_par_sets, pars_sd)
##' p_mat
##'
##' # Propose a new parameter set
##' p_mat <- pars$walk(p_mat, pars_sd)
##' p_mat
##'
##' # Information about parameters:
##' pars$names()
##' pars$summary()
##'
##' # Compute prior
##' pars$prior(p_mat)
##'
##' # Transform data for your model
##' pars$model(p_mat)
if2_parameters <- R6::R6Class(
  "if2_parameters",
  cloneable = FALSE,

  private = list(
    parameters = NULL,
    transform = NULL,
    integer = NULL,
    min = NULL,
    max = NULL
  ),

  public = list(
    ##' @description Create the if2_parameters object
    ##'
    ##' @param parameters A `list` of
    ##' [`mcstate::if2_parameter`] objects, each of which describe a
    ##' single parameter in your model. If `parameters` is named, then
    ##' these names must match the `$name` element of each parameter is
    ##' used (this is verified).
    ##'
    ##' @param transform An optional transformation function to apply
    ##' to your parameter vector immediately before passing it to the
    ##' model function. If not given, then [as.list] is
    ##' used, as dust models require this. However, if you need to
    ##' generate derived parameters from those being actively sampled
    ##' you can do arbitrary transformations here.
    initialize = function(parameters, transform = NULL) {
      parameters <- check_parameters(parameters, "if2_parameter")

      if (is.null(transform)) {
        transform <- as.list
      }
      assert_is(transform, "function")

      private$parameters <- parameters
      private$transform <- transform

      private$integer <- vlapply(private$parameters, "[[", "integer",
                                  USE.NAMES = FALSE)
      private$min <- vnapply(private$parameters, "[[", "min",
                             USE.NAMES = FALSE)
      private$max <- vnapply(private$parameters, "[[", "max",
                             USE.NAMES = FALSE)
    },

    ##' @description Return the initial parameter values as a named numeric
    ##' vector
    initial = function() {
      vnapply(private$parameters, "[[", "initial")
    },

    ##' @description Set up a parameter walk
    ##'
    ##' @param n_par_sets An integer number of parameter sets, which
    ##' defines the size of the population being peturbed.
    ##'
    ##' @param pars_sd A vector of standard deviations for the walk
    ##' of each parameter
    walk_initialise = function(n_par_sets, pars_sd) {
      n_par_sets <- assert_integer(n_par_sets)
      n_pars <- length(private$parameters)

      pars_mat <- matrix(self$initial(), n_pars, n_par_sets)
      rownames(pars_mat) <- names(self$names())
      self$walk(pars_mat, pars_sd)
    },

    ##' @description Propose a new parameter matrix given a current matrix
    ##' and walk standard deviation vector.
    ##'
    ##' @param pars A parameter matrix, from this function or
    ##' `$walk_initialise()`
    ##'
    ##' @param pars_sd A vector of standard deviations for the walk
    ##' of each parameter
    walk = function(pars, pars_sd) {
      stopifnot(length(pars_sd) == nrow(pars))
      n_par_sets <- ncol(pars)
      for (par_idx in seq_len(length(private$parameters))) {
        pars[par_idx, ] <-
          rnorm(n_par_sets, pars[par_idx, ], pars_sd[[par_idx]])
        if (private$integer[par_idx]) {
          pars[par_idx, ] <- vnapply(pars[par_idx, ], round)
        }
        pars[par_idx, pars[par_idx, ] < private$min[par_idx]] <-
          private$min[par_idx]
        pars[par_idx, pars[par_idx, ] > private$max[par_idx]] <-
          private$max[par_idx]
      }
      rownames(pars) <- self$names()
      pars
    },

    ##' @description Return the names of the parameters
    names = function() {
      names(private$parameters)
    },

    ##' @description Return a [`data.frame`] with information about
    ##' parameters (name, min, max, and integer).
    summary = function() {
      data_frame(name = self$names(),
                 min = private$min,
                 max = private$max,
                 discrete = private$integer,
                 integer = private$integer)
    },

    ##' @description Compute the prior for a parameter vector
    ##'
    ##' @param pars a parameter matrix from `$walk()`
    prior = function(pars) {
      n_pars <- length(private$parameters)
      stopifnot(nrow(pars) == n_pars)
      ret <- rep(0.0, ncol(pars))
      for (i in seq_len(n_pars)) {
        ret <- ret + private$parameters[[i]]$prior(pars[i, ])
      }
      ret
    },

    ##' @description Apply the model transformation function to a parameter
    ##' vector. Output is a list for lists, suitable for use with a dust
    ##' object with `pars_multi = TRUE`
    ##'
    ##' @param pars a parameter matrix from `$walk()`
    model = function(pars) {
      apply(pars, 2, private$transform)
    }
  ))
