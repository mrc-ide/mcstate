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
##' #Construct an object with two parameters:
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
    max = NULL,
    populations = NULL,
    varied = NULL
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
    ##' used, as dust models require this. However, if you need to
    ##' generate derived parameters from those being actively sampled
    ##' you can do arbitrary transformations here.
    ##'
    ##' @param populations Specifies the names of the different populations.
    ##' Ignored if no varied parameters are included.
    initialize = function(parameters, proposal, transform = NULL,
                          populations = NULL) {
      parameters <- check_parameters(parameters,
                      c("pmcmc_parameter", "pmcmc_varied_parameter"))

      # Split constructor for different proposal types. First asssumes array
      #  with 3 dimensions, second is matrix.
      varied <- vlapply(parameters, inherits, what = "pmcmc_varied_parameter")
      if (any(varied)) {
        # is.null(private$varied) quick check to see if all fixed
        private$varied <- vcapply(parameters[varied], "[[", "name")
        private$populations <- assert_character(populations)
        lapply(parameters, function(.x) {
          if (!(length(.x$initial) %in% c(1, length(populations)))) {
            stop(sprintf(
              "Expected length '1' or '%d' but got '%d' for parameter '%s'.",
              length(populations), length(.x$initial), .x$name
            ))
          }
        })

        # proposal
        assert_is(proposal, "array")
        dims <- dim(proposal)

        # assume that if matrix provided then same proposal for all populations
        if (is.na(dims[3])) {
          arr_proposal <- array(proposal,
            dim = c(dim(proposal), length(populations)))
          if (!is.null(dimnames(proposal))) {
            dimnames(arr_proposal) <- c(dimnames(proposal), list(populations))
          }
          proposal <- arr_proposal
          dims <- dim(proposal)
        }

        # catch wrong dimension sizes
        if (dims[3] != length(populations) ||
            !all(dims[1:2] == length(parameters))) {
          stop(sprintf(
            "Expected proposal array with dimensions %d x %d x %d.",
            length(parameters), length(parameters), length(populations)))
        }

        if (!is.null(dimnames(proposal))) {
          ## At this point we could reorder if that's useful and these
          ## are setequal.
          ok <- identical(rownames(proposal), names(parameters)) &&
            identical(colnames(proposal), names(parameters)) &&
            identical(dimnames(proposal)[3], populations)
          if (!ok) {
            stop("Expected dimension names of 'proposal' to match parameters
            and populations.")
          }
        }

         private$proposal <- apply(proposal, 3, rmvnorm_generator)
         private$min <- sapply(parameters, function(.x) {
           min <- .x[["min"]]
           if (length(min) == 1) {
             min <- rep(min, length(populations))
           }
           min
         })
         private$max <- sapply(parameters, function(.x) {
           max <- .x[["max"]]
           if (length(max) == 1) {
             max <- rep(max, length(populations))
           }
           max
         })
         rownames(private$min) <- populations
         rownames(private$max) <- populations

      } else {
        assert_is(proposal, "matrix")
        # proposal
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

         private$proposal <- rmvnorm_generator(proposal)
         private$min <- vnapply(parameters, "[[", "min",
                                USE.NAMES = FALSE)
         private$max <- vnapply(parameters, "[[", "max",
                                USE.NAMES = FALSE)
      }

      if (is.null(transform)) {
        transform <- as.list
      } else {
        assert_is(transform, "function")
      }

      private$parameters <- parameters
      private$proposal_kernel <- proposal
      private$transform <- transform

      private$discrete <- vlapply(private$parameters, "[[", "discrete",
                                  USE.NAMES = FALSE)
    },

    ##' @description Return the initial parameter values as a named numeric
    ##' vector or named matrix for multiple populations.
    # FIXME - COULD CONSIDER ALWAYS RETURNING MATRIX FOR CONSISTENCY
    initial = function() {
      if (length(private$populations)) {
        initial <- sapply(private$parameters, function(.x) {
           initial <- .x[["initial"]]
           if (length(initial) == 1) {
             initial <- rep(initial, length(private$populations))
           }
           initial
        })
        rownames(initial) <- private$populations
        initial
      } else {
        vnapply(private$parameters, "[[", "initial")
      }
    },

    ##' @description Return the names of the parameters
    names = function() {
      names(private$parameters)
    },

    ##' @description Return a `data.frame` with information about
    ##' parameters (name, min, max, discrete, fixed or varied).
    ##'
    ##' @param population For parameter sets including varied parameters,
    ##' `population` specifies which population to summarise. If `NULL` then
    ##'  returns summary of each population as a list. If no varied parameters
    ##'  included then `population` is ignored.
    summary = function(population = NULL) {
      if (is.null(private$varied)) {
        data_frame(
          name = self$names(),
          min = private$min,
          max = private$max,
          discrete = private$discrete
        )
      } else {
        if (is.null(population)) {
          # FIXME - not convinced if best representation, could consider array
          # or force 1 pop
          pops <- lapply(private$populations, function(.x) {
            data_frame(
              name = self$names(),
              min = private$min[rownames(private$min) %in% .x],
              max = private$max[rownames(private$max) %in% .x],
              discrete = private$discrete,
              type = ifelse(self$names() %in% private$varied, "varied", "fixed")
            )
          })
          names(pops) <- private$populations
          pops
        } else {
          if (!(population %in% private$populations)) {
            stop(sprintf("Expected 'population' in %s.",
                          str_collapse(private$populations)))
          }
          data_frame(
            name = self$names(),
            min = private$min[rownames(private$min) %in% population],
            max = private$max[rownames(private$max) %in% population],
            discrete = private$discrete,
            type = ifelse(self$names() %in% private$varied, "varied", "fixed")
          )
        }
      }
    },

    ##' @description Compute the prior(s) for a parameter vector/matrix
    ##'
    ##' @param theta a parameter vector in the same order as your
    ##' parameters were defined in (see `$names()` for that order.
    prior = function(theta) {
      if (is.null(private$populations)) {
        lp <- Map(function(p, value) p$prior(value), private$parameters, theta)
        sum(list_to_numeric(lp))
      } else {
        priors <- numeric(nrow(theta))
        names(priors) <- private$populations
        for (i in seq_along(priors)) {
          lp <- Map(function(p, value) {
            fprior <- p$prior
            if (is.function(fprior)) {
              fprior(value)
            } else {
              fprior[[i]](value)
            }
          }, private$parameters, theta[i, ])
          priors[[i]] <- sum(list_to_numeric(lp))
        }
        priors
      }
    },

    ##' @description Propose a new parameter vector given a current parameter
    ##' vector. This proposes a new parameter vector/matrix given your current
    ##' vector and the variance-covariance matrix of your proposal
    ##' kernel, discretises any discrete values, and reflects bounded
    ##' parameters until they lie within `min`:`max`. Returns matrix if any
    ##' varied parameters are included, otherwise vector.
    ##'
    ##' @param theta a parameter vector in the same order as your
    ##' parameters were defined in (see `$names()` for that order.
    ##'
    ##' @param scale an optional scaling factor to apply to the
    ##' proposal distribution. This may be useful in sampling starting
    ##' points. The parameter is equivalent to a multiplicative factor
    ##' applied to the variance covariance matrix.
    ##'
    ##' @param type specifies which type of parameters should be proposed,
    ##' either fixed parameters only ("fixed"), varied only ("varied"), or
    ##' both ("both") types. For 'fixed' and 'varied' then parameters of the
    ##' other type is returned equal to corresponding value in `theta`.
    # FIXME - AS WITH INITIAL COULD CONSIDER ALWAYS RETURNING MATRIX FOR
    #         CONSISTENCY
    propose = function(theta, scale = 1, type = c("fixed", "varied", "both")) {
      if (is.null(private$varied)) {
        theta <- private$proposal(theta, scale)
        theta[private$discrete] <- round(theta[private$discrete])
        reflect_proposal(theta, private$min, private$max)
      } else {
        proposals <- private$proposal
        type <- match.arg(type)
        mpropose <- matrix(nrow = nrow(theta), ncol = ncol(theta),
                           dimnames = dimnames(theta))
        for (i in seq_along(proposals)) {
          proposal <- proposals[[i]](theta[i, ], scale)
          proposal[private$discrete] <- round(proposal[private$discrete])
          proposal <- reflect_proposal(proposal, private$min[i, ],
                                       private$max[i, ])
          if (type != "both") {
              # if not requested type then revert to theta
              which <- self$summary(private$populations[[i]])$type != type
              proposal[which] <- theta[i, which]
          }
          mpropose[i, ] <- proposal
        }
        mpropose
      }
    },

    ##' @description Apply the model transformation function to a parameter
    ##' vector.
    ##'
    ##' @param theta a parameter vector in the same order as your
    ##' parameters were defined in (see `$names()` for that order.
    model = function(theta) {
      if (is.null(private$populations)) {
        private$transform(theta)
      } else {
        apply(theta, 1, private$transform)
      }
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

      base <- set_names(rep(NA_real_, length(private$parameters)),
                        names(private$parameters))
      base[idx_fixed] <- fixed
      base_transform <- private$transform

      if (is.null(private$populations)) {
        proposal <- private$proposal_kernel[idx_vary, idx_vary, drop = FALSE]
      } else {
        proposal <- array(
          apply(private$proposal_kernel, 3, "[", idx_vary,
                                idx_vary, drop = FALSE),
          dim = c(length(idx_vary), length(idx_vary),
                  length(private$populations)))
      }

      transform <- function(p) {
          base_transform(set_into(base, idx_vary, p))
      }

      pmcmc_parameters$new(private$parameters[idx_vary], proposal, transform,
        populations = private$populations)
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

# p <- pmcmc_parameters$new(
#   list(
#     pmcmc_parameter("beta", 0.2, min = 0, max = 1,
#                     prior = function(p) log(1e-10)),
#     pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
#                     prior = function(p) log(1e-10))),
#   proposal = proposal_kernel)

# pmcmc_parameters_shared$new(
#   list(
#     pmcmc_parameter("beta", 0.2, min = 0, max = 1,
#                     prior = function(p) log(1e-10)),
#     pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
#                     prior = function(p) log(1e-10))),
#   proposal = proposal_kernel,
#   shared = c("beta"),
#   group = regions)
## $propose gets argument type = c("individual", "shared", "both") and
## always returns a list
## $prior takes that list and returns a vector
## Assume 3 regions here
## individual:
## 1. propose individual parameters (varied)
##    [[beta + a, gamma], [beta + b, gamma], [beta + c, gamma]]
## 2. compute priors for all 3 (implies that the prior function
##    returns a vector)
## 3. compute likelihoods for each of the 3 regions, giving posteriors
## 4. accept or reject the updates individually
## shared:
## 1. propose shared parmeters
##    [[beta, gamma + x], [beta, gamma + x], [beta, gamma + x]]
## 2. compute sum log prior
## 3. compute likelihoods for each of the 3 regions, and sum
## 4. accept or reject collectively
## both:
## 1. propose all parameters
##    [[beta + a, gamma + x], [beta + b, gamma + x], [beta + c, gamma +x]]#
## 2. compute sum log prior
## 3. compute likelihoods for each of the 3 regions, and sum
## 4. accept or reject collectively