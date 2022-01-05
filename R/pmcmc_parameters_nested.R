##' @title pmcmc_parameters_nested
##'
##' @description Construct nested parameters for use with
##'   [pmcmc()]. This creates a utility object that is used
##'   internally to work with parameters that may be fixed and the same for all
##'   given populations, or varied and possibly-different between populations.
##'   Most users only need to construct this object, but see the examples for
##'   how it can be used.
##'
##' @export
##' @examples
##' # Construct an object with two varied parameters ('a' and 'b'),
##' # two fixed parameters ('c' and 'd') and two populations ('p1' and 'p2')
##' parameters <- list(mcstate::pmcmc_varied_parameter("a", c("p1", "p2"), 2),
##'                    mcstate::pmcmc_varied_parameter("b", c("p1", "p2"), 2),
##'                    mcstate::pmcmc_parameter("c", 3),
##'                    mcstate::pmcmc_parameter("d", 4))
##' proposal_fixed <- diag(2)
##' proposal_varied <- diag(2) + 1
##' pars <- mcstate::pmcmc_parameters_nested$new(parameters, proposal_varied,
##'                                              proposal_fixed)
##'
##' # Initial parameters
##' p <- pars$initial()
##' p
##'
##' # Propose a new parameter point
##' pars$propose(p, type = "both")
##' pars$propose(p, type = "fixed")
##' pars$propose(p, type = "varied")
##'
##' # Information about parameters:
##' pars$names()
##' pars$names("fixed")
##' pars$names("varied")
##' pars$summary()
##'
##' # Compute log prior probability, per population
##' pars$prior(p)
##'
##' # Transform data for your model
##' pars$model(p)
pmcmc_parameters_nested <- R6::R6Class(
  "pmcmc_parameters_nested",
  cloneable = FALSE,

  private = list(
    parameters = NULL,
    proposal_kernel = NULL,
    transform = NULL,
    inner = NULL,
    ## TODO: rename
    pops = NULL
  ),

  public = list(
    ##' @description Create the pmcmc_parameters object
    ##'
    ##' @param parameters A `list` of
    ##' [pmcmc_parameter] or [pmcmc_varied_parameter] objects, each of which
    ##' describe a single (possibly-varying) parameter in your model.
    ##' If `parameters` is named, then these names must match the `$name`
    ##' element of each parameter that is used (this is verified).
    ##'
    ##' @param proposal_varied,proposal_fixed Square proposal matrices
    ##' corresponding to the variance-covariance matrix of a multivariate
    ##' gaussian distribution used to generate new varied and fixed parameters
    ##' respectively.'. They must have the same number of
    ##' rows and columns as there are varied and fixed parameters respectively.
    ##' The names must correspond exactly to the names in
    ##' `parameters`. Because it corresponds to a variance-covariance
    ##' matrix it must be symmetric and positive definite.
    ##'
    ##' @param populations Specifies the names of the different populations
    ##' that the varying parameters change according to. Only required if no
    ##' [pmcmc_varied_parameter] objects are included in `parameters`.
    ##' Otherwise population names are taken from those objects.
    ##'
    ##' @param transform An optional transformation function to apply
    ##' to your parameter vector immediately before passing it to the
    ##' model function. If not given, then [as.list] is
    ##' used, as dust models require this. However, if you need to
    ##' generate derived parameters from those being actively sampled
    ##' you can do arbitrary transformations here.
    initialize = function(parameters, proposal_varied = NULL,
                          proposal_fixed = NULL, populations = NULL,
                          transform = NULL) {
      parameters <- ppn_validate_parameters(parameters)
      populations <- ppn_validate_populations(parameters, populations)
      proposal_kernel <- ppn_validate_proposals(parameters, proposal_varied,
                                                proposal_fixed)
      transform <- ppn_validate_transform(transform, populations)

      inner <- list()
      if (!is.null(parameters$fixed)) {
        inner$fixed <-
          pmcmc_parameters$new(parameters$fixed, proposal_kernel$fixed)
      }

      if (!is.null(parameters$varied)) {
        make_p <- function(p) {
          pmcmc_parameters$new(
            lapply(parameters$varied, "[[", p),
            array_drop(proposal_kernel$varied[, , p, drop = FALSE], 3))
        }
        inner$varied <- set_names(lapply(populations, make_p), populations)
      }

      private$parameters <- parameters
      private$pops <- populations
      private$proposal_kernel <- proposal_kernel
      private$inner <- inner
      private$transform <- transform
    },

    ##' @description Return the names of the parameters
    ##'
    ##' @param type One of "both" (the default, all parameters),
    ##'   "fixed" (parameters that are shared across populations) or
    ##'   "varied" (parameters that vary over populations).
    names = function(type = "both") {
      type <- match_value(type, c("both", "fixed", "varied"))
      names(private$parameters[[type]])
    },

    ##' @description Return the names of the populations
    populations = function() {
      private$pops
    },

    ##' @description Validate a parameter matrix.  This method
    ##' checks that your matrix has the expected size (rows according
    ##' to parameters, columns to populations) and if named that the
    ##' names are exactly what is expected.  It also verifies that the
    ##' fixed parameters are same across all populations.
    ##'
    ##' @param theta a parameter matrix
    validate = function(theta) {
      expected <- list(parameters = self$names(),
                       populations = self$populations())
      assert_dimensions(theta, lengths(expected, FALSE))
      theta <- assert_dimnames(theta, expected)

      i <- self$names("fixed")
      if (!all(theta[i, ] == theta[i, 1])) {
        stop("Fixed parameters are not everywhere fixed")
      }

      theta
    },

    ##' @description Return a `data.frame` with information about
    ##' parameters (name, min, max, discrete, type (fixed or varied)
    ##' and population)
    summary = function() {
      populations <- self$populations()
      fixed <- varied <- NULL
      if (length(self$names("fixed")) > 0) {
        fixed <- private$inner$fixed$summary()
        fixed <- cbind(
          fixed[rep(seq_len(nrow(fixed)), length(populations)), ],
          type = "fixed",
          population = rep(populations, each = nrow(fixed)),
          stringsAsFactors = FALSE)
      }

      if (length(self$names("varied")) > 0) {
        varied <- lapply(private$inner$varied, function(x) x$summary())
        varied <- cbind(
          do.call("rbind", varied),
          type = "varied",
          population = rep(populations, each = nrow(varied[[1]])),
          stringsAsFactors = FALSE)
      }

      ret <- rbind(fixed, varied)
      i <- order(match(ret$population, populations),
                 match(ret$name, self$names()))
      ret <- ret[i, ]
      rownames(ret) <- NULL
      ret
    },

    ##' @description Return the initial parameter values as a named matrix with
    ##' rows corresponding to parameters and columns to populations.
    initial = function() {
      pops <- self$populations()
      nms <- self$names()
      ret <- matrix(NA_real_, length(nms), length(pops),
                    dimnames = list(nms, pops))

      nms_varied <- self$names("varied")
      if (length(nms_varied) > 0) {
        ret[nms_varied, ] <-
          vapply(private$inner$varied, function(p) p$initial(),
                 numeric(length(nms_varied)))
      }

      nms_fixed <- self$names("fixed")
      if (length(nms_fixed) > 0) {
        ret[nms_fixed, ] <- private$inner$fixed$initial()
      }

      ret
    },

    ##' @description Compute the prior(s) for a parameter matrix. Returns a
    ##' named vector with names corresponding to populations.
    ##'
    ##' @param theta a parameter matrix with columns in the same order as
    ##' `$names()` and rows in the same order as `$populations()`.
    prior = function(theta) {
      theta <- self$validate(theta)
      pops <- self$populations()
      ret <- set_names(rep(0, length(pops)), pops)

      nms_fixed <- self$names("fixed")
      if (length(nms_fixed) > 0) {
        ret <- ret + private$inner$fixed$prior(theta[nms_fixed, 1])
      }

      nms_varied <- self$names("varied")
      if (length(nms_varied) > 0) {
        ret <- ret + vnapply(pops, function(x)
          private$inner$varied[[x]]$prior(theta[nms_varied, x]))
      }

      ret
    },

    ##' @description This proposes a new parameter matrix given your current
    ##' matrix and the variance-covariance matrices of the proposal
    ##' kernels, discretises any discrete values, and reflects bounded
    ##' parameters until they lie within `min`:`max`. Returns matrix with rows
    ##' corresponding to parameters and columns to populations (i.e.,
    ##' the same orientation as `theta`).
    ##'
    ##' @param theta a parameter matrix with rows in the same order as
    ##' `$names()` and columns in the same order as `$populations()`.
    ##'
    ##' @param type specifies which type of parameters should be proposed,
    ##' either fixed parameters only ("fixed"), varied only ("varied"), or
    ##' both ("both") types. For 'fixed' and 'varied', parameters of the
    ##' other type are left unchanged.
    ##'
    ##' @param scale an optional scaling factor to apply to the
    ##' proposal distribution. This may be useful in sampling starting
    ##' points. The parameter is equivalent to a multiplicative factor
    ##' applied to the variance covariance matrix.
    propose = function(theta, type, scale = 1) {
      theta <- self$validate(theta)
      type <- match_value(type, c("both", "varied", "fixed"))

      nms_fixed <- self$names("fixed")
      if (type %in% c("fixed", "both") && length(nms_fixed) > 0) {
        theta[nms_fixed, ] <-
          private$inner$fixed$propose(theta[nms_fixed, 1], scale)
      }

      nms_varied <- self$names("varied")
      if (type %in% c("varied", "both") && length(nms_varied) > 0) {
        theta[nms_varied, ] <-
          vapply(self$populations(), function(x)
            private$inner$varied[[x]]$propose(theta[nms_varied, x], scale),
            numeric(length(nms_varied)))
      }

      theta
    },

    ##' @description Apply the model transformation function to a parameter
    ##' matrix.
    ##'
    ##' @param theta a parameter matrix with rows in the same order as
    ##' `$names()` and columns in the same order as `$populations()`.
    model = function(theta) {
      nested_transform(self$validate(theta), private$transform)
    },

    ##' @description Set some parameters to fixed values. Use this to
    ##' reduce the dimensionality of your system.  Note that this function
    ##' has an unfortunate name collision - we use "fixed" and "varied"
    ##' parameters generally to refer to ones that are fixed across
    ##' populations or which vary among populations.  However, in the
    ##' context of this method "fixed" refers to parameters which will
    ##' be set to a single value and no longer used in inference.
    ##'
    ##' @param fixed a named vector of parameters to fix
    fix = function(fixed) {
      populations <- self$populations()
      fixed <- ppn_validate_fix(fixed, populations,
                                self$names("fixed"), self$names("varied"))

      i_fixed <- setdiff(self$names("fixed"), rownames(fixed))
      i_varied <- setdiff(self$names("varied"), rownames(fixed))
      i_keep <- setdiff(self$names(), rownames(fixed))

      parameters <- private$parameters$both[i_keep]

      if (length(i_fixed) > 0) {
        proposal_fixed <-
          private$proposal_kernel$fixed[i_fixed, i_fixed, drop = FALSE]
      } else {
        proposal_fixed <- NULL
      }

      if (length(i_varied) > 0) {
        proposal_varied <-
          private$proposal_kernel$varied[i_varied, i_varied, , drop = FALSE]
      } else {
        proposal_varied <- NULL
      }

      ## The transform management is hardest; need to provide a *list*
      ## of functions, one per population.
      make_transform <- function(pop) {
        base <- base[, pop, drop = TRUE]
        base_transform <- private$transform[[pop]]
        function(p) {
          base[i_keep] <- p
          base_transform(base)
        }
      }
      base <- matrix(NA_real_, length(self$names()), length(populations),
                     dimnames = list(self$names(), populations))
      base[rownames(fixed), ] <- fixed
      transform <- set_names(lapply(populations, make_transform), populations)

      pmcmc_parameters_nested$new(parameters,
                                  proposal_varied, proposal_fixed,
                                  populations, transform)
    }
  ))


## There are a great many things that have to align here, so we do
## them all at once so that it's not too boring.
ppn_validate_parameters <- function(parameters) {
  assert_is(parameters, "list")
  if (length(parameters) == 0) {
    stop("At least one parameter is required")
  }

  assert_list_of(parameters, c("pmcmc_parameter", "pmcmc_varied_parameter"))

  is_varied <- vlapply(parameters, inherits, "pmcmc_varied_parameter")

  fixed <- ppn_validate_parameters_fixed(parameters[!is_varied])
  varied <- ppn_validate_parameters_varied(parameters[is_varied])

  nms <- character(length(parameters))
  nms[!is_varied] <- names(fixed)
  nms[is_varied] <- names(varied)

  dups <- nms[duplicated(nms)]
  if (length(dups) > 0L) {
    stop("Duplicate parameter names: ",
         paste(squote(unique(dups)), collapse = ", "))
  }

  names(parameters) <- nms

  list(both = parameters, fixed = fixed, varied = varied)
}


ppn_validate_parameters_fixed <- function(parameters) {
  if (length(parameters) == 0) {
    return(NULL)
  }
  nms <- vcapply(parameters, "[[", "name", USE.NAMES = FALSE)
  if (!is.null(names(parameters)) && !identical(nms, names(parameters))) {
    stop("Fixed parameters are named, but the names do not match parameters")
  }
  names(parameters) <- nms
  parameters
}


ppn_validate_parameters_varied <- function(parameters) {
  if (length(parameters) == 0) {
    return(NULL)
  }

  ## only need to get the name of the first population as internally asserted
  nms <- vcapply(parameters, function(x) x[[1]]$name[1], USE.NAMES = FALSE)

  if (!is.null(names(parameters)) && !identical(nms, names(parameters))) {
    stop("Varied parameters are named, but the names do not match parameters")
  }
  names(parameters) <- nms

  populations <- lapply(parameters, function(x) names(x))
  ok <- vlapply(populations[-1], function(x) identical(x, populations[[1]]))
  if (!all(ok)) {
    stop("Populations and ordering of varied parameters must be identical")
  }

  parameters
}


ppn_validate_populations <- function(parameters, populations) {
  if (is.null(populations)) {
    if (is.null(parameters$varied)) {
      stop(paste("Either varied parameters must be included in",
                 "'parameters' or 'populations' must be non-NULL"))
    }
    populations <- names(parameters$varied[[1]])
  } else if (!is.null(parameters$varied)) {
    if (!identical(names(parameters$varied[[1]]), populations)) {
      stop("'population' does not match varied parameters")
    }
  }
  populations
}


ppn_validate_proposals <- function(parameters, proposal_varied,
                                   proposal_fixed) {
  if (is.null(parameters$varied)) {
    if (!is.null(proposal_varied)) {
      stop("'proposal_varied' supplied, but no varied parameters")
    }
  } else {
    if (is.null(proposal_varied)) {
      stop("'proposal_varied' not supplied for varied parameters")
    }
    populations <- names(parameters$varied[[1]])
    len_varied <- length(parameters$varied)
    len_pop <- length(populations)
    if (is.matrix(proposal_varied)) {
      assert_dimensions(proposal_varied, c(len_varied, len_varied))
      proposal_varied <- array(proposal_varied,
                               dim = c(len_varied, len_varied, len_pop))
    } else {
      assert_is(proposal_varied, "array")
      assert_dimensions(proposal_varied, c(len_varied, len_varied, len_pop))
    }
    proposal_varied <- assert_dimnames(
      proposal_varied,
      list(parameters = names(parameters$varied),
           parameters = names(parameters$varied),
           populations = populations))
  }

  if (is.null(parameters$fixed)) {
    if (!is.null(proposal_fixed)) {
      stop("'proposal_fixed' supplied, but no fixed parameters")
    }
  } else {
    if (is.null(proposal_fixed)) {
      stop("'proposal_fixed' not supplied for fixed parameters")
    }
    len_fixed <- length(parameters$fixed)
    assert_is(proposal_fixed, "matrix")
    assert_dimensions(proposal_fixed, c(len_fixed, len_fixed))
    proposal_fixed <- assert_dimnames(
      proposal_fixed,
      list(parameters = names(parameters$fixed),
           parameters = names(parameters$fixed)))
  }

  list(fixed = proposal_fixed, varied = proposal_varied)
}


ppn_validate_transform <- function(transform, populations) {
  if (is.null(transform)) {
    transform <- as.list
  }
  if (is.function(transform)) {
    transform <- set_names(rep(list(transform), length(populations)),
                           populations)
  } else {
    assert_list_of(transform, "function")
    if (!identical(names(transform), populations)) {
      stop("If 'transform' is a list, its names must be the populations")
    }
  }
  transform
}


ppn_validate_fix <- function(fixed, populations, names_fixed, names_varied) {
  assert_is(fixed, "matrix")
  names_all <- c(names_fixed, names_varied)
  if (!identical(colnames(fixed), populations)) {
    stop("colnames of 'fixed' must be identical to '$populations()'")
  }
  if (is.null(rownames(fixed))) {
    stop("'fixed' must have rownames (parameters)")
  }
  err <- setdiff(rownames(fixed), names_all)
  if (length(err) > 0) {
    stop("Fixed parameters not found in model: ",
         paste(squote(err), collapse = ", "))
  }
  if (anyDuplicated(rownames(fixed))) {
    stop("Duplicate fixed parameters")
  }
  if (nrow(fixed) == length(names_all)) {
    stop("Cannot fix all parameters")
  }

  i <- intersect(rownames(fixed), names_fixed)
  if (!all(fixed[i, ] == fixed[i, 1])) {
    ## This message is bad and I should feel bad.
    stop("Fixed fixed parameters are not everywhere fixed")
  }

  fixed
}


nested_transform <- function(theta, transform) {
  unname(lapply(seq_along(transform), function(i)
    transform[[i]](theta[, i])))
}
