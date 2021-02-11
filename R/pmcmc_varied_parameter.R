##' Describe a varying parameter for use within the nested pmcmc. Note that
##' the name is not set here, but will end up being naturally defined
##' when used with [`pmcmc_parameters_nested`], which collects
##' these together for use with [pmcmc()].
##'
##' @title Describe varying pmcmc parameter
##'
##' @param name Name for the parameter (a string)
##'
##' @param populations The name of the populations for which different values
##'   of the parameter are being estimated for, length `n_pop`.
##'
##' @param initial Initial value(s) for the parameter. Must be either length
##'   `n_pop` or `1`, in which case the same value is assumed for all
##'    populations.
##'
##' @param min Optional minimum value(s) for the parameter (otherwise
##'   `-Inf`). If given, then `initial` must be at least this
##'   value. Must be either length `n_pop` or `1`, in which case the same value
##'   is assumed for all populations.
##'
##' @param max Optional max value for the parameter (otherwise
##'   `Inf`). If given, then `initial` must be at most this
##'   value. Must be either length `n_pop` or `1`, in which case the same value
##'   is assumed for all populations.
##'
##' @param discrete Logical, indicating if this parameter is
##'   discrete. If `TRUE` then the parameter will be rounded
##'   after a new parameter is proposed. Must be either length `n_pop` or `1`,
##'   in which case the same value is assumed for all populations.
##'
##' @param prior A prior function (if not given an improper flat prior
##'   is used - be careful!). It must be a function that takes a
##'   single argument, being the value of this parameter. If given,
##'   then `prior(initial)` must evaluate to a finite value. Must be either
##'   length `n_pop` or `1`, in which case the same value is assumed for all
##'   populations.
##'
##' @export
##' @examples
##' pmcmc_varied_parameter(
##'   name = "size",
##'   populations = c("Europe", "America"),
##'   initial = c(100, 200),
##'   min = 0,
##'   max = Inf,
##'   discrete = TRUE,
##'   prior = list(dnorm, dexp)
##' )
pmcmc_varied_parameter <- function(name, populations, initial, min = -Inf,
                                   max = Inf, discrete = FALSE, prior = NULL) {

  assert_character(populations)
  len <- length(populations)

  initial <- recycle(initial, len)
  min <- recycle(min, len)
  max <- recycle(max, len)
  discrete <- recycle(discrete, len)

  if (is.null(prior)) {
    prior <- rep(list(function(p) 0), len)
  } else {
    if (is.function(prior)) {
      prior <- rep(list(prior), len)
    } else {
      assert_is(prior, "list")
      prior <- recycle(prior, len)
    }
  }

  params <- Map(pmcmc_parameter, name, initial, min, max, discrete, prior)
  names(params) <- populations
  class(params) <- "pmcmc_varied_parameter"

  params
}
