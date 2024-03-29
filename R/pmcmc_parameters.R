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
##' @param discrete  Deprecated; use `integer` instead.
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
##' @param mean Optionally, an estimate of the mean of the
##'   parameter. If not given, then we assume that `initial` is a
##'   reasonable estimate.  This is used only in adaptive mcmc.
##'
##' @export
##' @examples
##' pmcmc_parameter("a", 0.1)
pmcmc_parameter <- function(name, initial, min = -Inf, max = Inf,
                            discrete, integer = FALSE, prior = NULL,
                            mean = NULL) {
  if (!missing(discrete)) {
    .Deprecated("integer", old = "discrete")
    integer <- discrete
  }
  assert_scalar_character(name)
  assert_scalar_logical(integer)
  assert_in_range(initial, min, max)

  if (is.null(mean)) {
    mean <- initial
  }
  assert_in_range(mean, min, max)

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
              integer = integer, prior = prior, mean = mean)
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
##' @section Parameter transformations:
##'
##' Unless you have a very simple model, it is highly unlikely that
##'   the parameters that you are interested in performing inference
##'   on are the same as the parameters that you might need to
##'   initialise your model.
##'
##' Due to the nature of mcmc and other inference algorithms, the
##'   general assumption is that the inference parameters will be a
##'   simple vector of real values; here each of the `parameters`
##'   elements corresponds to one of these. The proposal matrix maps
##'   one vector to another via a simple multivariate-gaussian kernel.
##'
##' On the other hand, dust models can take a named list of arbitrary
##'   data as their input parameters (see
##'   [dust::dust_generator]). These might include:
##'
##' * things that are not parameters at all from the perspective of
##'   the inference - for example some quantity that you might vary
##'   depending on the region/species/etc you're running the model for
##'   but that you are not fitting.
##' * non-scalar quantities that are directly derived from some
##'   parameters that you are fitting.  As an example of this, in
##'   [sircovid](https://mrc-ide.github.io), a transmission model of
##'   COVID, we take a number of "contact rates" which apply at
##'   different points in time, and generate from this an interpolated
##'   series of contact rates per time step (a very long
##'   vector). Other users have needed to generate equilibrium
##'   solutions to parts of their model and used these at
##'   initialisation.
##' * arbitrary complex inputs to the model, for example weather data,
##'   demographic matrices, population contact rate matrices
##'   etc. These are all "parameters" from the perspective of a dust
##'   model but not at all from the perspective of the inference
##'   process.
##'
##' To allow for this in a flexible way, mcstate allows a "transform"
##'   function, the `transform` argument to the constructor. This
##'   function maps a named numeric vector of inference parameters to
##'   whatever you need for your dust model.  The default value for
##'   this function is [as.list] which just converts the named vector
##'   to a named list, which works well in the example cases here.
##'
##' When providing a transformation function, you may want to provide
##'   a "closure" rather than a top-level function. This way you can
##'   bind additional data into your function.  For example, suppose
##'   that you want to use some demographic matrix `m` in your model,
##'   and perform inference on parameters `a` and `b` you might write
##'
##' ```
##' make_transform <- function(m) {
##'   function(theta) {
##'     c(list(m = m), as.list(theta))
##'   }
##' }
##' ```
##'
##' and pass this into `mcstate::pmcmc_parameters$new`, providing
##'   parameter definitions only for `a` and `b`.  See the examples
##'   for full working of this.
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
##'
##' # Above we describe a nontrivial transformation function using a closure
##' make_transform <- function(m) {
##'   function(theta) {
##'     c(list(m = m), as.list(theta))
##'   }
##' }
##'
##' # Suppose this is our demographic matrix (note here that the name
##' # need not match that used in the transform)
##' demographic_matrix <- diag(4)
##'
##' # Construct the parameters as above, but this time passing in the
##' # function that make_transform returns
##' pars <- mcstate::pmcmc_parameters$new(
##'   list(mcstate::pmcmc_parameter("a", 0.1, min = 0, max = 1,
##'                                 prior = function(a) log(a)),
##'        mcstate::pmcmc_parameter("b", 0, prior = dnorm)),
##'   matrix(c(1, 0.5, 0.5, 2), 2, 2),
##'   make_transform(demographic_matrix))
##'
##' # Now, as above we start from a position in terms of a and b only:
##' pars$initial()
##'
##' # But when prepared for the model, our matrix will be set up
##' pars$model(pars$initial())
pmcmc_parameters <- R6::R6Class(
  "pmcmc_parameters",
  cloneable = FALSE,

  private = list(
    parameters = NULL,
    proposal = NULL,
    proposal_kernel = NULL,
    transform = NULL,
    integer = NULL,
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

      assert_is(proposal, "matrix")
      if (!all(dim(proposal) == length(parameters))) {
        stop(sprintf(
          "Expected a square proposal matrix with %d rows and columns",
          length(parameters)))
      }
      if (!is.null(dimnames(proposal))) {
        ## At this point we could reorder if that's useful and these
        ## are setequal.
        ok <- identical(rownames(proposal), names(parameters)) &&
          identical(colnames(proposal), names(parameters))
        if (!ok) {
          stop("Expected dimension names of 'proposal' to match parameters")
        }
      }

      private$parameters <- parameters
      private$proposal_kernel <- proposal
      private$proposal <- rmvnorm_generator(proposal)
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

    ##' @description Return the estimate of the mean of the parameters,
    ##'   as set when created (this is not updated by any fitting!)
    mean = function() {
      vnapply(private$parameters, "[[", "mean")
    },

    ##' @description Return the variance-covariance matrix used for the
    ##'   proposal.
    vcv = function() {
      private$proposal_kernel
    },

    ##' @description Return the names of the parameters
    names = function() {
      names(private$parameters)
    },

    ##' @description Return a `data.frame` with information about
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
    ##' @param theta a parameter vector in the same order as your
    ##' parameters were defined in (see `$names()` for that order.
    prior = function(theta) {
      lp <- Map(function(p, value) p$prior(value), private$parameters, theta)
      sum(list_to_numeric(lp))
    },

    ##' @description Propose a new parameter vector given a current parameter
    ##' vector. This proposes a new parameter vector given your current
    ##' vector and the variance-covariance matrix of your proposal
    ##' kernel, rounds any integer values, and reflects bounded
    ##' parameters until they lie within `min`:`max`.
    ##'
    ##' @param theta a parameter vector in the same order as your
    ##' parameters were defined in (see `$names()` for that order.
    ##'
    ##' @param scale an optional scaling factor to apply to the
    ##' proposal distribution. This may be useful in sampling starting
    ##' points. The parameter is equivalent to a multiplicative factor
    ##' applied to the variance covariance matrix.
    ##'
    ##' @param vcv A variance covariance matrix of the correct size,
    ##' overriding the proposal matrix built into the parameters object.
    ##' This will be slightly less efficient but allow a different proposal
    ##' matrix to be used (e.g., during an adaptive MCMC)
    propose = function(theta, scale = 1, vcv = NULL) {
      if (is.null(vcv)) {
        theta_new <- private$proposal(theta, scale)
      } else {
        theta_new <- rmvnorm_generator(vcv, check = FALSE)(theta)
      }
      theta_new[private$integer] <- round(theta_new[private$integer])
      reflect_proposal(theta_new, private$min, private$max)
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
