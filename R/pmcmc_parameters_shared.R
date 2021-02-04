##' @title pmcmc_parameters_shared
##'
##' @description Construct possibly-shared parameters for use with
##'   [pmcmc()]. This creates a utility object that is used
##'   internally to work with possibly-shared parameters. Most users only
##'   need to construct this object, but see the examples for how it can be
##'   used.
##'
##' @export
##' @examples
##' #Construct an object with four parameters, two varied, two fixed:
##' pars <- mcstate::pmcmc_parameters$new(
##'   list(mcstate::pmcmc_parameter("a", 0.1, min = 0, max = 1,
##'                                 prior = function(a) log(a)),
##'        mcstate::pmcmc_parameter("b", 0, prior = dnorm),
##'        mcstate::pmcmc_parameter("c", 0.5, min = 0, max = 1,
##'                                 prior = dunif),
##'        mcstate::pmcmc_parameter("d", 1)),
##'   proposal = matrix(c(1, 0.5, 0.5, 2), 2, 2),
##'   varied = c("a", "c"),
##'   populations = c("Europe", "America")
##' )
pmcmc_parameters_shared <- R6::R6Class(
  "pmcmc_parameters_shared",
  cloneable = FALSE,
  inherit = pmcmc_parameters,

  private = list(
    parameters = NULL,
    proposal = NULL,
    proposal_kernel = NULL,
    transform = NULL,
    discrete = NULL,
    min = NULL,
    max = NULL,
    varied = NULL,
    fixed = NULL,
    populations = NULL
  ),

  public = list(
    ##' @description Create the pmcmc_parameters_shared object
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
    ##'
    ##' @param varied Optional character vector to specify which parameters
    ##' can vary between populations. If not given then assumed all parameters
    ##' are fixed.
    ##'
    ##' @param populations Specifies the names of the different populations.
    ##' Ignored if `varied` is NULL.
    ##
    ## TODO:
    ##  - Need to think about if transform can stay as-is or needs to be updated.
    ##  - Starting with simple case of same proposal for each region, likely will
    ##      need to overload constructor if this is extended to list/array.
    initialize = function(parameters, proposal, transform = NULL,
    varied = NULL, populations = NULL) {
    super$initialize(parameters, proposal, transform)

    par_names <- self$names()

    # Check 'varied' in parameter names then add fixed as set difference
    # 'populations' ignored if varied is NULL (not required to do this)
    if (!is.null(varied)) {
        if (!all(varied %in% par_names)) {
            stop(sprintf(
                "Not all 'varied' in parameter names (%s)",
                str_collapse(par_names)
            ))
        } else {
            private$varied <- varied
            private$fixed <- setdiff(par_names, varied)
            private$populations <- assert_character(populations)
        }
    }

    invisible(self)

    },

    ##' @description Return the initial parameter values as a named numeric
    ##' vector
    initial = function() {
      vnapply(private$parameters, "[[", "initial")
    },

    ##' @description Return a `data.frame` with information about
    ##' parameters (name, min, max, and discrete).
    summary = function() {
      data_frame(
          name = self$names(),
          min = private$min,
          max = private$max,
          discrete = private$discrete,
          type = ifelse(self$names() %in% private$varied, "varied", "fixed")
      )
    },

    ##' @description Compute the priors for a parameter list
    ##'
    ##' @param theta a parameter list in the same order as your
    ##' parameters were defined in (see `$names()` for that order.
    ##'
    ##' @param type specifies which type of parameters should be proposed,
    ##' either fixed parameters only ("fixed"), varied only ("varied"), or
    ##' both ("both") types. For 'fixed' and 'varied' then parameters of the
    ##' other type is returned equal to corresponding value in `theta`.
    prior = function(theta, type = c("fixed", "varied", "both")) {
        sapply(theta, function(.x) {
            lp <- Map(function(p, value) p$prior(value), private$parameters, .x)
            sum(list_to_numeric(lp))
        })
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
    ##'
    ##' @param type specifies which type of parameters should be proposed,
    ##' either fixed parameters only ("fixed"), varied only ("varied"), or
    ##' both ("both") types. For 'fixed' and 'varied' then parameters of the
    ##' other type is returned equal to corresponding value in `theta`.
    ##
    ## TODO - Whilst currently only single proposal accepted, written as if
    ##   a list.
    propose = function(theta, scale = 1, type = c("fixed", "varied", "both")) {
      pmat <- private$proposal
      type <- match.arg(type)
      # convert to list if only one matrix provided for all regions
      # FIXME - Though this may be done in construction so could delete below
      if (!inherits(pmat, "list")) {
          pmat <- rep(list(pmat), length(private$populations))
      }

      lapply(pmat, function(.x) {
        proposal <- .x(theta, scale)
        proposal[private$discrete] <- round(proposal[private$discrete])
        proposal <- reflect_proposal(proposal, private$min, private$max)
        # if proposal not requested for both then revert proposal to theta for
        # the 'other' type
        if (type != "both") {
            which <- self$summary()$type != type
            proposal[which] <- theta[which]
        }
        proposal
      })
    }
  ))


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