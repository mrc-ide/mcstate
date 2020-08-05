pmcmc_parameter <- function(initial, min = -Inf, max = Inf, discrete = FALSE,
                            prior = NULL) {
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
        stop("Failed to evalute your prior function on initial value: ",
             e$message))
    if (!is.finite(value)) {
      stop("Your prior function returned an non-finite value on initial value")
    }
  }

  ret <- list(initial = initial, min = min, max = max, discrete = discrete,
              prior = prior %||% function(p) 0)
  class(ret) <- "pmcmc_parameter"
  ret
}


pmcmc_parameters <- R6::R6Class(
  "pmcmc_parameters",

  private = list(
    parameters = NULL,
    proposal = NULL,
    transform = NULL,
    discrete = NULL,
    min = NULL,
    max = NULL
  ),

  public = list(
    initialize = function(parameters, proposal, transform = NULL) {
      ok <- vlapply(parameters, inherits, "pmcmc_parameter")
      if (!all(ok)) {
        stop("Expected all elements of '...' to be 'pmcmc_parameter' objects")
      }

      expect_named(parameters)

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
          stop("Expected dimension names of 'proposal' to match parmeters")
        }
      }

      private$parameters <- parameters
      private$proposal <- rmvnorm_generator(proposal)
      private$transform <- transform

      private$discrete <- vlapply(private$parameters, "[[", "discrete")
      private$min <- vnapply(private$parameters, "[[", "min")
      private$max <- vnapply(private$parameters, "[[", "max")
    },

    initial = function() {
      vnapply(private$parameters, "[[", "initial")
    },

    prior = function(theta) {
      lp <- Map(function(p, value) p$prior(value), private$parameters, theta)
      sum(list_to_numeric(lp))
    },

    propose = function(theta) {
      theta <- private$proposal(theta)
      theta[private$discrete] <- round(theta[private$discrete])
      reflect_proposal(theta, private$min, private$max)
    },

    model = function(theta) {
      private$transform(theta)
    }
  ))


## create function to reflect proposal boundaries at pars_min and pars_max
## this ensures the proposal is symetrical and we can simplify the MH step
reflect_proposal <- function(x, x_min, x_max) {
  ## TODO: Cope with infinite or missing ranges
  x_range <- x_max - x_min
  abs((x + x_range - x_min) %% (2 * x_range) - x_range) + x_min
}
