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
