##' Describe a single parameter for use within the SMC^2. Note that
##' the name is not set here, but will end up being naturally defined
##' when used with [`smc2_parameters`], which collects
##' these together for use with [smc2()].
##'
##' @title Describe single pmcmc parameter
##'
##' @param name Name for the parameter (a string)
##'
##' @param sample A sampling function; it must take a single argument
##'   representing the number of sampled to be returned. Typically
##'   this will be a `r` probability function corresponding to the
##'   sampling version of your prior (e.g., you might use `runif` and
##'   `dunif` for `sample` and `prior`). If you provide `min`, `max`
##'   or `discrete` you *must* ensure that your function returns
##'   values that satisfy these constraints, as this is not (yet)
##'   checked.
##'
##' @param prior A prior function. It must be a function that takes a
##'   single argument, being the value of this parameter.
##'
##' @param min Optional minimum value for the parameter (otherwise
##'   `-Inf`). If given, then `initial` must be at least this
##'   value.
##'
##' @param max Optional max value for the parameter (otherwise
##'   `Inf`). If given, then `initial` must be at most this
##'   value.
##'
##' @param discrete Logical, indicating if this parameter is
##'   discrete. If `TRUE` then the parameter will be rounded
##'   after a new parameter is proposed.
##'
##' @export
##' @examples
##' mcstate::smc2_parameter("a",
##'                         function(n) rnorm(n),
##'                         function(x) dnorm(n, log = TRUE))
smc2_parameter <- function(name, sample, prior,
                           min = -Inf, max = Inf, discrete = FALSE) {
  assert_scalar_character(name)
  assert_scalar_logical(discrete)
  assert_is(sample, "function")
  assert_is(prior, "function")

  if (max <= min) {
    stop(sprintf("'max' must be > 'min' (%s)", min))
  }

  ## TODO: Would Raphael's distribution classes help us out here for
  ## the priors? There's a natural coupling that is awkward here.
  ##
  ## https://github.com/alan-turing-institute/distr6
  ##
  ## Overhead looks like a factor of 500-1000x over the raw functions,
  ## but it's not going to be a massive timesink, and it would clean
  ## up the interface nicely.
  ret <- list(name = name, sample = sample, prior = prior,
              min = min, max = max, discrete = discrete)
  class(ret) <- "smc2_parameter"
  ret
}


##' @title pmcmc_parameters
##'
##' @description Construct parameters for use with
##'   [smc()]. This creates a utility object that is used
##'   internally to work with parameters. Most users only need to
##'   construct this object, but see the examples for how it can be
##'   used.
##'
##' @export
smc2_parameters <- R6::R6Class(
  "smc2_parameters",
  cloneable = FALSE,

  private = list(
    parameters = NULL,
    transform = NULL,
    discrete = NULL,
    min = NULL,
    max = NULL,

    constrain_parameters = function(theta) {
      theta[, private$discrete] <- round(theta[, private$discrete])
      min <- rep(private$min, each = nrow(theta))
      max <- rep(private$max, each = nrow(theta))
      theta[] <- reflect_proposal(theta, min, max)
      theta
    }
  ),

  public = list(
    ##' @description Create the smc2_parameters object
    ##'
    ##' @param parameters A `list` of
    ##' [smc2_parameter] objects, each of which describe a
    ##' single parameter in your model. If `parameters` is named, then
    ##' these names must match the `$name` element of each parameter is
    ##' used (this is verified).
    ##'
    ##' @param transform An optional transformation function to apply
    ##' to your parameter vector immediately before passing it to the
    ##' model function. If not given, then [as.list] is
    ##' used, as dust models require this. However, if t you need to
    ##' generate derived parameters from those being actively sampled
    ##' you can do arbitrary transformations here.
    initialize = function(parameters, transform = NULL) {
      parameters <- check_parameters(parameters, "smc2_parameter")

      if (is.null(transform)) {
        transform <- as.list
      }
      assert_is(transform, "function")

      private$parameters <- parameters
      private$transform <- transform

      private$discrete <- vlapply(private$parameters, "[[", "discrete",
                                  USE.NAMES = FALSE)
      private$min <- vnapply(private$parameters, "[[", "min",
                             USE.NAMES = FALSE)
      private$max <- vnapply(private$parameters, "[[", "max",
                             USE.NAMES = FALSE)
    },

    ##' @description Create `n` independent random parameter vectors (as a
    ##'   matrix with `n` rows)
    ##'
    ##' @param n Number of replicate parameter sets to draw
    sample = function(n) {
      ret <- vapply(private$parameters, function(p) p$sample(n), numeric(n))
      private$constrain_parameters(ret)
    },

    ##' @description Return the names of the parameters
    names = function() {
      names(private$parameters)
    },

    ##' @description Return a `data.frame` with information about
    ##' parameters (name, min, max, and discrete).
    summary = function() {
      data_frame(name = self$names(),
                 min = private$min,
                 max = private$max,
                 discrete = private$discrete)
    },

    ##' @description Compute the prior for a parameter vector
    ##'
    ##' @param theta a parameter vector in the same order as your
    ##' parameters were defined in (see `$names()` for that order.
    prior = function(theta) {
      np <- length(private$parameters)
      stopifnot(ncol(theta) == np)
      ret <- 0.0
      for (i in seq_len(np)) {
        ret <- ret + private$parameters[[i]]$prior(theta[, i])
      }
      ret
    },

    ##' @description Propose a new parameter vector given a current parameter
    ##' vector and variance covariance matrix. After proposal, this discretises
    ##' any discrete values, and reflects bounded parameters until they lie
    ##' within `min`:`max`.
    ##'
    ##' @param theta a parameter vector in the same order as your
    ##' parameters were defined in (see `$names()` for that order).
    ##'
    ##' @param vcv the variance covariance matrix for the proposal; must
    ##' be square and have a number of rows and columns equal to the
    ##' number of parameters, in the same order as `theta`.
    propose = function(theta, vcv) {
      np <- length(private$parameters)
      stopifnot(ncol(theta) == np,
                nrow(vcv) == np,
                ncol(vcv) == np)
      private$constrain_parameters(t(apply(theta, 1, rmvnorm_generator(vcv))))
    },

    ##' @description Apply the model transformation function to a parameter
    ##' vector.
    ##'
    ##' @param theta a parameter vector in the same order as your
    ##' parameters were defined in (see `$names()` for that order.
    model = function(theta) {
      lapply(seq_len(nrow(theta)), function(i)
             private$transform(theta[i, ]))
    }
  ))
