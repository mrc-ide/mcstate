##' Describe a single parameter for use within the pmcmc. Note that
##' the name is not set here, but will end up being naturally defined
##' when used with [`pmcmc_parameters`], which collects
##' these together for use with [pmcmc()].
##'
##' @title Describe single pmcmc parameter
##'
##' @param name Name for the parameter (a string)
##'
##' @param initial Initial value for the parameter
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
##' @param prior A prior function (if not given an improper flat prior
##'   is used - be careful!). It must be a function that takes a
##'   single argument, being the value of this parameter. If given,
##'   then `prior(initial)` must evaluate to a finite value.
##'
##' @export
##' @examples
##' pmcmc_parameter("a", 0.1)
pmcmc_parameter <- function(name, initial, min = -Inf, max = Inf,
                            discrete = FALSE, prior = NULL) {
  assert_scalar_character(name)
  assert_scalar_logical(discrete)

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
              discrete = discrete, prior = prior)
  class(ret) <- "pmcmc_parameter"
  ret
}


##' @title pmcmc_parameters
##'
##' @description Construct parameters for use with
##'   [pmcmc()]. This creates a utility object that is used
##'   internally to work with parameters. Most users only need to
##'   construct this object, but see the examples for how it can be
##'   used.
##'
##' @export
##' @examples
##' # Construct an object with two parameters:
##' pars <- mcstate::pmcmc_parameters$new(
##'   list(mcstate::pmcmc_parameter("a", 0.1, min = 0, max = 1,
##'                                 prior = function(a) log(a)),
##'        mcstate::pmcmc_parameter("b", 0, prior = dnorm)),
##'   matrix(c(1, 0.5, 0.5, 2), 2, 2))
##'
##' # Initial parameters
##' p <- pars$initial()
##' p
##'
##' # Propose a new parameter point
##' pars$propose(p)
##'
##' # Information about parameters:
##' pars$names()
##' pars$summary()
##'
##' # Compute prior
##' pars$prior(p)
##'
##' # Transform data for your model
##' pars$model(p)
pmcmc_parameters <- R6::R6Class(
  "pmcmc_parameters",
  cloneable = FALSE,

  private = list(
    parameters = NULL,
    proposal = NULL,
    proposal_kernel = NULL,
    transform = NULL,
    discrete = NULL,
    min = NULL,
    max = NULL
  ),

  public = list(
    ##' @description Create the pmcmc_parameters object
    ##'
    ##' @param parameters A `list` of
    ##' [pmcmc_parameter] objects, each of which describe a
    ##' single parameter in your model. If `parameters` is named, then
    ##' these names must match the `$name` element of each parameter is
    ##' used (this is verified).
    ##'
    ##' @param proposal A square proposal distribution corresponding to the
    ##' variance-covariance matrix of a multivariate gaussian distribution
    ##' used to generate new parameters. It must have the same number of
    ##' rows and columns as there are elements in `parameters`, and if
    ##' named the names must correspond exactly to the names in
    ##' `parameters`. Because it corresponds to a variance-covariance
    ##' matrix it must be symmetric and positive definite.
    ##'
    ##' @param transform An optional transformation function to apply
    ##' to your parameter vector immediately before passing it to the
    ##' model function. If not given, then [as.list] is
    ##' used, as dust models require this. However, if t you need to
    ##' generate derived parameters from those being actively sampled
    ##' you can do arbitrary transformations here.
    initialize = function(parameters, proposal, transform = NULL) {
      parameters <- check_parameters(parameters, "pmcmc_parameter")

      if (is.null(transform)) {
        transform <- as.list
      }
      assert_is(transform, "function")

      private$parameters <- parameters
      private$transform <- transform

      self$update_proposal(proposal)

      private$discrete <- vlapply(private$parameters, "[[", "discrete",
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
      lp <- Map(function(p, value) p$prior(value), private$parameters, theta)
      sum(list_to_numeric(lp))
    },

    ##' @description Propose a new parameter vector given a current parameter
    ##' vector. This proposes a new parameter vector given your current
    ##' vector and the variance-covariance matrix of your proposal
    ##' kernel, discretises any discrete values, and reflects bounded
    ##' parameters until they lie within `min`:`max`.
    ##'
    ##' @param theta a parameter vector in the same order as your
    ##' parameters were defined in (see `$names()` for that order.
    ##'
    ##' @param scale an optional scaling factor to apply to the
    ##' proposal distribution. This may be useful in sampling starting
    ##' points. The parameter is equivalent to a multiplicative factor
    ##' applied to the variance covariance matrix.
    propose = function(theta, scale = 1) {
      theta <- private$proposal(theta, scale)
      theta[private$discrete] <- round(theta[private$discrete])
      reflect_proposal(theta, private$min, private$max)
    },

    ##' Update the proposal kernel.  Use this to adapt the proposal kernel
    ##' originally set into the object mid-run.
    ##'
    ##' @param proposal A square proposal distribution corresponding to the
    ##' variance-covariance matrix of a multivariate gaussian distribution
    ##' used to generate new parameters (the same restrictions apply as
    ##' for `proposal` passed to the constructor).
    update_proposal = function(proposal) {
      assert_is(proposal, "matrix")
      n <- length(private$parameters)
      if (!all(dim(proposal) == c(n, n))) {
        stop(sprintf(
          "Expected a square proposal matrix with %d rows and columns", n))
      }
      if (!is.null(dimnames(proposal))) {
        expected <- list(names(private$parameters), names(private$parameters))
        if (!identical(unname(dimnames(proposal)), expected)) {
          stop("Expected dimension names of 'proposal' to match parameters")
        }
      }
      private$proposal_kernel <- proposal
      private$proposal <- rmvnorm_generator(proposal)
    },

    ##' @description Apply the model transformation function to a parameter
    ##' vector.
    ##'
    ##' @param theta a parameter vector in the same order as your
    ##' parameters were defined in (see `$names()` for that order.
    model = function(theta) {
      private$transform(theta)
    },

    ##' @description Set some parameters to fixed values. Use this to
    ##' reduce the dimensionality of your system.
    ##'
    ##' @param fixed a named vector of parameters to fix
    fix = function(fixed) {
      assert_named(fixed, TRUE)
      idx_fixed <- match(names(fixed), names(private$parameters))
      if (any(is.na(idx_fixed))) {
        stop("Fixed parameters not found in model: ",
             paste(squote(names(fixed)[is.na(idx_fixed)]), collapse = ", "))
      }
      if (length(idx_fixed) == length(private$parameters)) {
        stop("Cannot fix all parameters")
      }
      idx_vary <- setdiff(seq_along(private$parameters), idx_fixed)
      proposal <- private$proposal_kernel[idx_vary, idx_vary, drop = FALSE]

      base <- set_names(rep(NA_real_, length(private$parameters)),
                        names(private$parameters))
      base[idx_fixed] <- fixed
      base_transform <- private$transform
      transform <- function(p) {
        base_transform(set_into(base, idx_vary, p))
      }
      pmcmc_parameters$new(private$parameters[idx_vary], proposal, transform)
    }
  ))


## create function to reflect proposal boundaries at pars_min and pars_max
## this ensures the proposal is symetrical and we can simplify the MH step
reflect_proposal <- function(x, x_min, x_max) {
  i <- x < x_min | x > x_max
  if (any(i)) {
    i_both <- i & is.finite(x_min) & is.finite(x_max)
    i_min <- i & is.finite(x_min) & !is.finite(x_max)
    i_max <- i & !is.finite(x_min) & is.finite(x_max)
    x[i_both] <- reflect_proposal_both(x[i_both], x_min[i_both], x_max[i_both])
    x[i_min] <- reflect_proposal_one(x[i_min], x_min[i_min])
    x[i_max] <- reflect_proposal_one(x[i_max], x_max[i_max])
  }
  x
}


reflect_proposal_both <- function(x, x_min, x_max) {
  x_r <- x_max - x_min
  abs((x + x_r - x_min) %% (2 * x_r) - x_r) + x_min
}


reflect_proposal_one <- function(x, x_bound) {
  2 * x_bound - x
}


check_parameters <- function(parameters, type) {
  assert_is(parameters, "list")
  if (length(parameters) == 0) {
    stop("At least one parameter is required")
  }
  ok <- vlapply(parameters, inherits, type)
  if (!all(ok)) {
    stop(sprintf("Expected all elements of '...' to be '%s' objects", type))
  }
  nms <- vcapply(parameters, "[[", "name", USE.NAMES = FALSE)
  dups <- nms[duplicated(nms)]
  if (length(dups) > 0L) {
    stop("Duplicate parameter names: ",
         paste(squote(unique(dups)), collapse = ", "))
  }
  if (!is.null(names(parameters)) && !identical(nms, names(parameters))) {
    stop("'parameters' is named, but the names do not match parameters")
  }
  names(parameters) <- nms
  parameters
}
