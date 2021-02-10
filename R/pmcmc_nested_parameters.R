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
pmcmc_nested_parameters <- R6::R6Class(
  "pmcmc_nested_parameters",
  cloneable = FALSE,

  private = list(
    fixed_parameters = NULL,
    varied_parameters = NULL,
    proposal_varied = NULL,
    proposal_fixed = NULL,
    proposal_kernel_varied = NULL,
    proposal_kernel_fixed = NULL,
    transform = NULL,
    param_names = NULL
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
    ##' @param populations Specifies the names of the different populations
    ##' that the parameters vary according to. Assumes that varied parameters
    ##' are estimated based on representative samples drawn from the respective
    ##' populations. Ignored if no varied parameters are included.
    initialize = function(parameters, proposal_fixed, proposal_varied,
      transform = NULL) {

      parameters <- clean_parameters(parameters)
      kernels <- clean_proposals(proposal_fixed, proposal_varied,
        parameters$pop, parameters$names)

      # avoids constant looping later
      private$param_names <- parameters$names

      # create pmcmc_parameters
      private$fixed_parameters <-
        pmcmc_parameters$new(parameters$fixed, kernels$fixed)

      private$varied_parameters <- vector("list", length(parameters$pop))
      names(private$varied_parameters) <- parameters$pop
      for (i in seq_along(parameters$pop)) {
        kernel <- kernels$varied[, , i]
        # catch for dropping on numeric
        if (!inherits(kernel, "matrix")) {
          kernel <- as.matrix(kernel)
        }
        private$varied_parameters[[i]] <-
         pmcmc_parameters$new(
           lapply(parameters$varied, "[[", parameters$pop[[i]]), kernel)
      }

      proposals <- make_proposals(kernels$fixed, kernels$varied,
        parameters$names)

      private$proposal_kernel_fixed <- proposals$kernel_fixed
      private$proposal_kernel_varied <- proposals$kernel_varied
      private$proposal_fixed <- proposals$proposal_fixed
      private$proposal_varied <- proposals$proposal_varied
      names(private$proposal_kernel_varied) <- parameters$pop
      names(private$proposal_varied) <- parameters$pop

      # transform function
      if (is.null(transform)) {
        private$transform <- as.list
      } else {
        private$transform <- assert_is(transform, "function")
      }

      invisible(self)
    },

    ##' @description Return the names of the parameters
    names = function() {
      unlist(private$param_names, use.names = FALSE)
    },

    ##' @description Return the names of the populations
    populations = function() {
      names(private$varied_parameters)
    },

        ##' @description Return a `data.frame` with information about
    ##' parameters (name, min, max, discrete, fixed or varied).
    ##'
    ##' @param population For parameter sets including varied parameters,
    ##' `population` specifies which population to summarise. If `NULL` then
    ##'  returns summary of each population as a list.
    summary = function(population = NULL) {
      summarise_population <- function(population) {
        if (!(population %in% self$populations())) {
          stop(sprintf("Expected 'population' in %s.",
            str_collapse(self$populations())))
        }
        type <- rep("fixed", length(self$names()))
        type[self$names() %in% private$param_names$varied] <- "varied"
        min <- c(r6_private(private$fixed_parameters)$min,
          r6_private(private$varied_parameters[[population]])$min)
        max <- c(r6_private(private$fixed_parameters)$max,
          r6_private(private$varied_parameters[[population]])$max)
        discrete <- c(r6_private(private$fixed_parameters)$discrete,
          r6_private(private$varied_parameters[[population]])$discrete)

        data_frame(
          name = self$names(),
          min = min,
          max = max,
          discrete = discrete,
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

    ##' @description Return the initial parameter values as a named numeric
    ##' vector or named matrix for multiple populations.
    initial = function() {
      fixed <- private$fixed_parameters$initial()
      varied <- vapply(private$varied_parameters, function(x) x$initial(),
        numeric(length(self$populations())))

      t(rbind(
        matrix(fixed, 2, 2, dimnames = list(names(fixed), self$populations())),
        varied))
    },

    ##' @description Compute the prior(s) for a parameter vector/matrix
    ##'
    ##' @param theta a parameter matrix with columns in the same order as
    ##' `$names()` and rows in the same order as `$populations()`.
    prior = function(theta) {
      priors <- numeric(nrow(theta))
      names(priors) <- self$populations()

      lp_fix <- private$fixed_parameters$prior(
        theta[1, private$param_names$fixed])
      lp_vary <- vnapply(self$populations(),
        function(x) private$varied_parameters[[x]]$prior(
          theta[x, private$param_names$varied]))

      lp_fix + lp_vary
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
    propose = function(theta, scale = 1, type = c("fixed", "varied", "both")) {
      type <- match.arg(type)
      pfix <- function(theta, scale) {
        private$proposal_fixed(theta, scale)
      }
      pvar <- function(theta, scale) {
        var <- vapply(seq(nrow(theta)), function(i)
          private$proposal_varied[[rownames(theta)[[i]]]](theta[i, ], scale),
          numeric(length(self$names())))
        colnames(var) <- self$populations()
        t(var)
      }
      if (type == "fixed") {
        pfix(theta, scale)
      } else if (type == "varied") {
        pvar(theta, scale)
      } else {
        fix <- pfix(theta, scale)
        var <- pvar(theta, scale)
        fix[, private$param_names$varied] <- var[, private$param_names$varied]
        fix
      }
    },

    ##' @description Apply the model transformation function to a parameter
    ##' vector.
    ##'
    ##' @param theta a parameter vector in the same order as your
    ##' parameters were defined in (see `$names()` for that order.
    model = function(theta) {
      apply(theta, 1, private$transform)
    },

    ##' @description Set some parameters to fixed values. Use this to
    ##' reduce the dimensionality of your system.
    ##'
    ##' @param fixed a named vector of parameters to fix
    # TODO
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

      if (is.null(self$populations())) {
        proposal <- private$proposal_kernel[idx_vary, idx_vary, drop = FALSE]
      } else {
        proposal <- array(
          apply(private$proposal_kernel, 3, "[", idx_vary,
                                idx_vary, drop = FALSE),
          dim = c(length(idx_vary), length(idx_vary),
                  length(self$populations())))
      }

      transform <- function(p) {
          base_transform(set_into(base, idx_vary, p))
      }

      pmcmc_parameters$new(private$parameters[idx_vary], proposal, transform,
        populations = self$populations())
    }
  ))


clean_parameters <- function(parameters) {
  assert_is(parameters, "list")
  if (length(parameters) == 0) {
    stop("At least one parameter is required")
  }

  type <- c("pmcmc_parameter", "pmcmc_varied_parameter")
  ok <- vlapply(parameters, inherits, type)
  if (!all(ok)) {
    stop(sprintf("Expected all elements of '...' to be '%s' objects", type))
  }

  varied <- vlapply(parameters, inherits, "pmcmc_varied_parameter")
  if (!any(varied)) {
    stop("No varied parameters, use `pmcmc_parameters` instead.")
  } else if (all(varied)) {
    fixed_parameters <- NULL
  } else {
    fixed_parameters <- parameters[!varied]
  }
  varied_parameters <- parameters[varied]

  # check fixed and varied names separately
  if (!is.null(fixed_parameters)) {
    fix_nms <- vcapply(fixed_parameters, "[[", "name", USE.NAMES = FALSE)
  } else {
    fix_nms <- c()
  }

  # only need to get the name of the first population as internally asserted
  var_nms <- vcapply(varied_parameters, function(x) x[[1]]$name[1],
    USE.NAMES = FALSE)
  nms <- c(fix_nms, var_nms)
  dups <- nms[duplicated(nms)]
  if (length(dups) > 0L) {
    stop("Duplicate parameter names: ",
         paste(squote(unique(dups)), collapse = ", "))
  }

  if (!is.null(fixed_parameters)) {
    if (!is.null(names(fixed_parameters)) &&
      !identical(fix_nms, names(fixed_parameters))) {
      stop("Fixed parameters are named, but the names do not match
      parameters.")
    }
    names(fixed_parameters) <- fix_nms
  }

  if (!is.null(names(varied_parameters)) &&
    !identical(var_nms, names(varied_parameters))) {
    stop("Varied parameters are named, but the names do not match parameters.")
  }
  names(varied_parameters) <- var_nms

  populations <- names(varied_parameters[[1]])
  ok <- vlapply(varied_parameters,
    function(x) identical(names(x), populations))
  if (!all(ok)) {
    stop("Populations and ordering of varied parameters should be identical.")
  }

  list(fixed = fixed_parameters, varied = varied_parameters, pop = populations,
    names = list(fixed = fix_nms, varied = var_nms))
}

clean_proposals <- function(proposal_fixed, proposal_varied, populations,
  params) {
  # minimal checks for proposals as these will be passed to constructors
  assert_is(proposal_fixed, "matrix")
  assert_is(proposal_varied, "array")

  if (is.na(dim(proposal_varied)[3])) {
     proposal_varied <- array(proposal_varied,
       dim = c(dim(proposal_varied), length(populations)))
  }

  if (!identical(dim(proposal_varied),
     c(length(params$varied), length(params$varied), length(populations)))) {
      stop(sprintf(
        "Expected proposal array with dimensions %d x %d x %d.",
        length(params$varied), length(params$varied), length(populations)))
  }

  if (!is.null(dimnames(proposal_varied)[[3]])) {
    if (!identical(dimnames(proposal_varied)[[3]], populations)) {
      stop("Expected 3rd dimension names of 'proposal' to match populations.")
    }
  }
  dimnames(proposal_varied)[3] <- list(populations)

  if (!is.null(rownames(proposal_varied)) ||
    !is.null(colnames(proposal_varied))) {
    ok <- identical(rownames(proposal_varied), params$varied) &&
            identical(colnames(proposal_varied), params$varied)
    if (!ok) {
      stop("Expected row and column names of 'proposal_varied' to match
      varied parameters.")
    }
  }
  dimnames(proposal_varied)[1:2] <- list(params$varied, params$varied)

  if (!is.null(rownames(proposal_fixed)) ||
    !is.null(colnames(proposal_fixed))) {
    ok <- identical(rownames(proposal_fixed), params$fixed) &&
            identical(colnames(proposal_fixed), params$fixed)
    if (!ok) {
      stop("Expected row and column names of 'proposal_fixed' to match
      fixed parameters.")
    }
  }
  dimnames(proposal_fixed)[1:2] <- list(params$fixed, params$fixed)


  list(fixed = proposal_fixed, varied = proposal_varied)
}

make_proposals <- function(kernel_fixed, kernel_varied, param_names) {
  kernel_varied_list <- lapply(seq_len(dim(kernel_varied)[[3]]),
    function(i) {
      kernel <- kernel_varied[, , i]
      # catch for drop on numeric
      if (!inherits(kernel, "matrix")) {
        kernel <- matrix(kernel,
          dimnames = rep(list(rownames(kernel_varied)), 2))
      }
      expand(kernel, unlist(param_names))
    })

  proposal_varied <- lapply(kernel_varied_list, rmvnorm_generator)

  proposal_fixed <- make_proposal_fixed(kernel_fixed)

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
