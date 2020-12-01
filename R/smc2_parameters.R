smc2_parameter <- function(name, sample, prior,
                           min = -Inf, max = Inf, discrete = FALSE) {
  assert_scalar_character(name)
  assert_scalar_logical(discrete)
  assert_is(prior, "function")

  if (max <= min) {
    stop(sprintf("'max' must be > 'min' (%s)", min))
  }

  ## TODO: We assume that the sample function automatically takes care
  ## of both the min/max and the discretisation here?
  ret <- list(name = name, sample = sample, prior = prior,
              min = min, max = max, discrete = discrete)
  class(ret) <- "smc2_parameter"
  ret
}


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

    sample = function(n) {
      ret <- vapply(private$parameters, function(p) p$sample(n), numeric(n))
      private$constrain_parameters(ret)
    },

    names = function() {
      names(private$parameters)
    },

    summary = function() {
      data_frame(name = self$names(),
                 min = private$min,
                 max = private$max,
                 discrete = private$discrete)
    },

    prior = function(theta) {
      np <- length(private$parameters)
      stopifnot(ncol(theta) == np)
      ret <- 0.0
      for (i in seq_len(np)) {
        ret <- ret + private$parameters[[i]]$prior(theta[, i])
      }
      ret
    },

    propose = function(theta, vcv) {
      np <- length(private$parameters)
      stopifnot(ncol(theta) == np,
                nrow(vcv) == np,
                ncol(vcv) == np)
      ## This could be made more efficient and generatie but I do not think that
      ## it's totally needed
      private$constrain_parameters(t(apply(theta, 1, rmvnorm_generator(vcv))))
    },

    model = function(theta) {
      lapply(seq_len(nrow(theta)), function(i)
             private$transform(theta[i, ]))
    }
  ))
