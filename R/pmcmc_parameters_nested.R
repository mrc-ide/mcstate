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
##' #Construct an object with two varied parameters, two fixed parameters,
##' #and two populations:
##' parameters <- list(pmcmc_varied_parameter("a", c("p1", "p2"), 2),
##'                    pmcmc_varied_parameter("b", c("p1", "p2"), 2),
##'                    pmcmc_parameter("c", 3),
##'                    pmcmc_parameter("d", 4))
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
pmcmc_parameters_nested <- R6::R6Class(
  "pmcmc_parameters_nested",
  cloneable = FALSE,

  private = list(
    fixed_parameters = list(),
    varied_parameters = list(),
    proposal_varied = list(),
    proposal_fixed = list(),
    proposal_kernel_varied = list(),
    proposal_kernel_fixed = list(),
    transform = NULL,
    param_names = NULL,
    pops = NULL,
    parameters = list()
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

      # original format required for `fix`
      private$parameters <- parameters

      parameters <- clean_parameters(parameters,
                                     proposal_varied, proposal_fixed)

      if (is.null(populations)) {
        if (is.null(proposal_varied)) {
          stop("Either varied parameters must be included in 'parameters' or
          'populations' must be non-NULL")
        } else {
          populations <- parameters$pop
        }
      }

      private$pops <- populations

      kernels <- clean_proposals(proposal_varied, proposal_fixed,
                                 populations, parameters$names)
      proposals <- make_proposals(kernels$varied, kernels$fixed,
                                  parameters$names)

      private$param_names <- parameters$names

      # create pmcmc_parameters
      if (!is.null(kernels$fixed)) {
        private$fixed_parameters <-
          pmcmc_parameters$new(parameters$fixed, kernels$fixed)

        private$proposal_kernel_fixed <- proposals$kernel_fixed
        private$proposal_fixed <- proposals$proposal_fixed
      }

      if (!is.null(kernels$varied)) {
        private$varied_parameters <- vector("list", length(populations))
        names(private$varied_parameters) <- populations
        for (i in seq_along(populations)) {
          kernel <- kernels$varied[, , i]
          # catch for dropping on numeric
          if (!inherits(kernel, "matrix")) {
            kernel <- as.matrix(kernel)
          }
          private$varied_parameters[[i]] <-
            pmcmc_parameters$new(
              lapply(parameters$varied, "[[", populations[[i]]), kernel)
        }

        private$proposal_kernel_varied <- proposals$kernel_varied
        private$proposal_varied <- proposals$proposal_varied
        names(private$proposal_kernel_varied) <- populations
        names(private$proposal_varied) <- populations
      }

      # transform function
      if (is.null(transform)) {
        private$transform <- function(x) apply(x, 1, as.list)
      } else {
        private$transform <- assert_is(transform, "function")
      }

      invisible(self)
    },

    ##' @description Return the names of the parameters
    names = function() {
      private$param_names$original
    },

    ##' @description Return the names of the populations
    populations = function() {
      private$pops
    },

    ##' @description Return a `data.frame` with information about
    ##' parameters (name, min, max, discrete, fixed or varied).
    ##'
    ##' @param population Specifies which population to summarise. If `NULL`
    ##' then returns a summary of each population as a list.
    summary = function(population = NULL) {
      summarise_population <- function(population) {
        if (!(population %in% self$populations())) {
          stop(sprintf("Expected 'population' in %s",
                       str_collapse(self$populations())))
        }
        type <- rep("fixed", length(self$names()))
        type[self$names() %in% private$param_names$varied] <- "varied"

        data_frame(
          name = self$names(),
          min = get_inner_params("min", population, private),
          max = get_inner_params("max", population, private),
          discrete = get_inner_params("discrete", population, private),
          type = type
        )
      }

      if (is.null(population)) {
        sum <- lapply(self$populations(), summarise_population)
        names(sum) <- self$populations()
        sum
      } else {
        summarise_population(population)
      }
    },

    ##' @description Return the initial parameter values as a named matrix with
    ##' rows corresponding to populations and columns are parameters.
    initial = function() {

      varied_params <- private$param_names$varied
      fixed_params <- private$param_names$fixed
      pops <- self$populations()

      if (length(varied_params) > 0) {
        varied <- matrix(vapply(private$varied_parameters,
                                function(x) x$initial(),
                                numeric(length(varied_params))),
                         ncol = length(pops),
                         nrow = length(private$param_names$varied),
                         dimnames = list(varied_params, pops))
      } else {
        varied <- NULL
      }

      if (length(fixed_params) > 0) {
        fixed <- private$fixed_parameters$initial()
        fixed <- matrix(fixed, length(fixed_params), length(pops),
                        dimnames = list(fixed_params, pops))
      } else {
        fixed <- NULL
      }

      init <- t(rbind(fixed, varied))
      init[, match(self$names(), colnames(init)), drop = FALSE]
    },

    ##' @description Compute the prior(s) for a parameter matrix. Returns a
    ##' named vector with names corresponding to populations.
    ##'
    ##' @param theta a parameter matrix with columns in the same order as
    ##' `$names()` and rows in the same order as `$populations()`.
    prior = function(theta) {
      pops <- self$populations()
      priors <- numeric(nrow(theta))
      names(priors) <- names(priors)

      if (length(private$fixed_parameters) > 0) {
        ## as theta is same in all populations for fixed, 1 chosen arbitrarily
        lp_fix <- private$fixed_parameters$prior(
          theta[1, private$param_names$fixed])
        # rep isn't necessary because of built-in sum vectorisation but it
        #   catches the no-varied case (same as names)
        lp_fix <- set_names(rep(lp_fix, length(pops)), pops)
      } else {
        lp_fix <- 0
      }

      if (length(private$varied_parameters) > 0) {
        lp_vary <- vnapply(pops,
                           function(x) private$varied_parameters[[x]]$prior(
                             theta[x, private$param_names$varied]))
      } else {
        lp_vary <- 0
      }

      lp_fix + lp_vary
    },

    ##' @description This proposes a new parameter matrix given your current
    ##' matrix and the variance-covariance matrices of the proposal
    ##' kernels, discretises any discrete values, and reflects bounded
    ##' parameters until they lie within `min`:`max`. Returns matrix with rows
    ##' corresponding to populations and columns to parameters.
    ##'
    ##' @param theta a parameter matrix with columns in the same order as
    ##' `$names()` and rows in the same order as `$populations()`.
    ##'
    ##' @param scale an optional scaling factor to apply to the
    ##' proposal distribution. This may be useful in sampling starting
    ##' points. The parameter is equivalent to a multiplicative factor
    ##' applied to the variance covariance matrix.
    ##'
    ##' @param type specifies which type of parameters should be proposed,
    ##' either fixed parameters only ("fixed"), varied only ("varied"), or
    ##' both ("both") types. For 'fixed' and 'varied', parameters of the
    ##' other type are left unchanged.
    propose = function(theta, scale = 1, type = c("fixed", "varied", "both")) {
      type <- match.arg(type)
      fix_params <- private$param_names$fixed
      var_params <- private$param_names$varied

      pfix <- function(theta, scale) {
        if (length(private$proposal_fixed) > 0) {
          private$proposal_fixed(theta, scale)
        }
      }
      pvar <- function(theta, scale) {
        if (length(private$proposal_varied) > 0) {
          var <- vapply(seq(nrow(theta)), function(i)
            private$proposal_varied[[rownames(theta)[[i]]]](theta[i, ], scale),
            numeric(length(self$names())))
          # catch one varied edge case
          if (!inherits(var, "matrix")) {
            var <- matrix(var, nrow = length(var_params),
                          dimnames = list(var_params, NULL))
          }

          colnames(var) <- self$populations()
          t(var)
        }
      }

      if (type == "fixed") {
        ret <- pfix(theta, scale)
      } else if (type == "varied") {
        ret <- pvar(theta, scale)
      } else {
        fix <- pfix(theta, scale)
        var <- pvar(theta, scale)
        if (is.null(fix)) {
          ret <- var
        } else if (is.null(var)) {
          ret <- fix
        } else {
          fix[, var_params] <- var[, var_params]
          ret <- fix
        }
      }

      if (!is.null(ret)) {
        discrete_params <- character(0)
        if (length(fix_params) > 0) {
          which <- r6_private(private$fixed_parameters)$discrete
          discrete_params <- c(discrete_params, fix_params[which])
        }
        if (length(var_params) > 0) {
          # as discrete is same across populations arbitrarily select first
          which <- r6_private(private$varied_parameters[[1]])$discrete
          discrete_params <- c(discrete_params, var_params[which])
        }

        ret[, discrete_params] <- round(ret[, discrete_params])

        for (i in seq(nrow(ret))) {
          min <- vnapply(private$parameters, function(x)
            if (is.null(x$min)) x[[i]]$min else x$min)
          max <- vnapply(private$parameters, function(x)
            if (is.null(x$max)) x[[i]]$max else x$max)
          ret[i, ] <- reflect_proposal(ret[i, ], min, max)
        }
      }

      ret
    },

    ##' @description Apply the model transformation function to a parameter
    ##' matrix.
    ##'
    ##' @param theta a parameter matrix with columns in the same order as
    ##' `$names()` and rows in the same order as `$populations()`.
    model = function(theta) {
      private$transform(theta)
    },

    ##' @description Set some parameters to fixed values. Use this to
    ##' reduce the dimensionality of your system.
    ##'
    ##' @param fixed a named vector of parameters to fix
    fix = function(fixed) {
      assert_is(fixed, "matrix")
      if (is.null(rownames(fixed)) || is.null(colnames(fixed))) {
        stop("'fixed' should have rownames (populations) and colnames
             (parameters)")
      }

      param_names <- private$param_names$original
      pops <- self$populations()

      if (any(is.na(match(rownames(fixed), pops)))) {
        stop("Rownames of 'fixed' should be identical to '$populations'")
      }

      idx_fixed <- match(colnames(fixed), param_names)
      if (any(is.na(idx_fixed))) {
        stop("Fixed parameters not found in model: ",
             paste(squote(colnames(fixed)[is.na(idx_fixed)]), collapse = ", "))
      }
      if (length(idx_fixed) == length(param_names)) {
        stop("Cannot fix all parameters")
      }
      idx_vary <- setdiff(seq_along(param_names), idx_fixed)

      base <- matrix(NA_real_, nrow = length(pops), ncol = length(param_names),
                     dimnames = list(pops, param_names))
      base[, idx_fixed] <- fixed
      base_transform <- private$transform

      idfix_fixed <- match(colnames(fixed), private$param_names$fixed)
      idfix_varied <- setdiff(seq_along(private$param_names$fixed),
                              idfix_fixed)
      idvar_fixed <- match(colnames(fixed), private$param_names$varied)
      idvar_varied <- setdiff(seq_along(private$param_names$varied),
                              idvar_fixed)
      if (length(idfix_varied) > 0) {
        proposal_fixed <- private$proposal_kernel_fixed[idfix_varied,
                                                        idfix_varied,
                                                        drop = FALSE]
      } else {
        proposal_fixed <- NULL
      }

      if (length(idvar_varied) > 0) {
        idvars <- match(private$param_names$varied,
                        private$param_names$original)
        proposal_varied <- simplify2array(private$proposal_kernel_varied)
        proposal_varied <- proposal_varied[idvars, idvars, , drop = FALSE]
        proposal_varied <- proposal_varied[idvar_varied, idvar_varied, ,
                                           drop = FALSE]
      } else {
        proposal_varied <- NULL
      }

      transform <- function(p) {
        base[, idx_vary] <- p
        base_transform(base)
      }

      pmcmc_parameters_nested$new(private$parameters[idx_vary],
                                  proposal_varied, proposal_fixed,
                                  self$populations(), transform)
    }
  ))


clean_parameters <- function(parameters, proposal_varied, proposal_fixed) {
  assert_is(parameters, "list")
  if (length(parameters) == 0) {
    stop("At least one parameter is required")
  }

  type <- c("pmcmc_parameter", "pmcmc_varied_parameter")
  ok <- vlapply(parameters, inherits, type)
  if (!all(ok)) {
    stop(sprintf(
      "Expected all elements of '...' to be '%s' objects",
      str_collapse(type)
    ))
  }

  varied <- vlapply(parameters, inherits, "pmcmc_varied_parameter")
  fixed_parameters <- parameters[!varied]
  varied_parameters <- parameters[varied]
  # convert empty list to NULL
  if (!length(fixed_parameters)) {
    fixed_parameters <- NULL
  } else {
    if (is.null(proposal_fixed)) {
      stop("'proposal_fixed' not supplied for fixed parameters")
    }
  }
  if (!length(varied_parameters)) {
    varied_parameters <- NULL
  } else {
    if (is.null(proposal_varied)) {
      stop("'proposal_varied' not supplied for fixed parameters")
    }
  }

  # check fixed and varied names separately
  fix_nms <- vcapply(fixed_parameters, "[[", "name", USE.NAMES = FALSE)

  # only need to get the name of the first population as internally asserted
  var_nms <- vcapply(varied_parameters, function(x) x[[1]]$name[1],
                     USE.NAMES = FALSE)
  nms <- c(fix_nms, var_nms)
  dups <- nms[duplicated(nms)]
  if (length(dups) > 0L) {
    stop("Duplicate parameter names: ",
         paste(squote(unique(dups)), collapse = ", "))
  }

  # store both original ordering and separated fixed/varied for reference
  orig_nms <- character(length(nms))
  orig_nms[!varied] <- fix_nms
  orig_nms[varied] <- var_nms

  if (!is.null(fixed_parameters)) {
    if (!is.null(names(fixed_parameters)) &&
        !identical(fix_nms, names(fixed_parameters))) {
      stop("Fixed parameters are named, but the names do not match
      parameters")
    }
    names(fixed_parameters) <- fix_nms
  }

  if (!is.null(varied_parameters)) {
    if (!is.null(names(varied_parameters)) &&
        !identical(var_nms, names(varied_parameters))) {
      stop("Varied parameters are named, but the names do not match
            parameters.")
    }
    names(varied_parameters) <- var_nms

    populations <- names(varied_parameters[[1]])
    ok <- vlapply(varied_parameters,
                  function(x) identical(names(x), populations))
    if (!all(ok)) {
      stop("Populations and ordering of varied parameters should be
            identical")
    }
  } else {
    populations <- NULL
  }

  list(fixed = fixed_parameters, varied = varied_parameters, pop = populations,
       names = list(fixed = fix_nms, varied = var_nms, original = orig_nms))
}


clean_proposals <- function(proposal_varied, proposal_fixed, populations,
                            params) {

  if (!is.null(proposal_varied)) {
    len_varied <- length(params$varied)
    len_pop <- length(populations)
    dim_names_var <- dimnames(proposal_varied)

    assert_is(proposal_varied, "array")

    if (is.na(dim(proposal_varied)[3])) {
      proposal_varied <- array(proposal_varied,
                               dim = c(dim(proposal_varied), len_pop))
    }

    if (!identical(dim(proposal_varied), c(len_varied, len_varied, len_pop))) {
      stop(sprintf(
        "Expected proposal array with dimensions %d x %d x %d",
        len_varied, len_varied, len_pop))
    }

    if (length(dim_names_var) == 3) {
      if (!identical(dim_names_var[[3]], populations)) {
        stop("Expected 3rd dimension names of 'proposal' to match
        populations")
      }
    }
    dimnames(proposal_varied)[3] <- list(populations)

    if (!is.null(dim_names_var[[1]]) || !is.null(dim_names_var[[2]])) {
      ok <- identical(dim_names_var[[1]], params$varied) &&
        identical(dim_names_var[[2]], params$varied)
      if (!ok) {
        stop("Expected row and column names of 'proposal_varied' to match
        varied parameters")
      }
    }
    dimnames(proposal_varied)[1:2] <- list(params$varied, params$varied)
  }

  # minimal checks for fixed proposals as these will be passed to constructors
  if (!is.null(proposal_fixed)) {
    dim_names_fix <- dimnames(proposal_fixed)
    assert_is(proposal_fixed, "matrix")
    if (!is.null(dim_names_fix[[1]]) || !is.null(dim_names_fix[[2]])) {
      ok <- identical(dim_names_fix[[1]], params$fixed) &&
        identical(dim_names_fix[[2]], params$fixed)
      if (!ok) {
        stop("Expected row and column names of 'proposal_fixed' to match
        fixed parameters")
      }
    }
    dimnames(proposal_fixed)[1:2] <- list(params$fixed, params$fixed)
  }

  list(fixed = proposal_fixed, varied = proposal_varied)
}


make_proposals <- function(kernel_varied, kernel_fixed, params) {
  if (!is.null(kernel_varied)) {
    kernel_varied_list <- lapply(seq_len(dim(kernel_varied)[[3]]), function(i)
      expand(kernel_varied[, , i, drop = FALSE], params$original))
    names(kernel_varied_list) <- dimnames(kernel_varied)[[3]]
    proposal_varied <- lapply(kernel_varied_list, rmvnorm_generator)
  } else {
    proposal_varied <- NULL
    kernel_varied_list <- NULL
  }

  if (!is.null(kernel_fixed)) {
    proposal_fixed <- make_proposal_fixed(kernel_fixed)
  } else {
    proposal_fixed <- NULL
  }

  list(kernel_varied = kernel_varied_list,
       kernel_fixed = kernel_fixed,
       proposal_varied = proposal_varied,
       proposal_fixed = proposal_fixed)

}


expand <- function(m, complete) {
  n <- length(complete)
  ret <- matrix(0, ncol = n, nrow = n, dimnames = list(complete, complete))
  ret[rownames(m), colnames(m)] <- m
  ret
}


make_proposal_fixed <- function(kernel) {
  rmv <- rmvnorm_generator(kernel)
  nms <- rownames(kernel)
  function(mean, scale = 1.0) {
    mean[, nms] <- rep(rmv(mean[1, nms], scale = scale), each = nrow(mean))
    mean
  }
}


get_inner_params <- function(param, population, private) {
  if (length(private$fixed_parameters) > 0) {
    fixed <- r6_private(private$fixed_parameters)[[param]]
    names(fixed) <- private$param_names$fixed
  } else {
    fixed <- NULL
  }

  if (length(private$varied_parameters) > 0) {
    varied <- r6_private(private$varied_parameters[[population]])[[param]]
    names(varied) <- private$param_names$varied
  } else {
    varied <- NULL
  }

  unname(c(fixed, varied)[private$param_names$original])
}
