pmcmc_varied_parameter <- function(name, populations, initial, min = -Inf,
                                  max = Inf, discrete = FALSE, prior = NULL) {

  assert_character(populations)
  lng <- length(populations)

  if (length(initial) == 1) {
    initial <- rep(initial, lng)
  } else if (length(initial) != lng) {
    stop(sprintf("Length of 'initial' must be '1' or %d.", lng))
  }

  if (length(min) == 1) {
    min <- rep(min, lng)
  } else if (length(min) != lng) {
    stop(sprintf("Length of 'min' must be '1' or %d.", lng))
  }

  if (length(max) == 1) {
    max <- rep(max, lng)
  } else if (length(max) != lng) {
    stop(sprintf("Length of 'max' must be '1' or %d.", lng))
  }

  if (length(discrete) == 1) {
    discrete <- rep(discrete, lng)
  } else if (length(discrete) != lng) {
    stop(sprintf("Length of 'discrete' must be '1' or %d.", lng))
  }


  if (is.null(prior)) {
    prior <- rep(list(function(p) 0), lng)
  } else {
    if (is.function(prior)) {
      prior <- rep(list(prior), lng)
    } else {
      assert_is(prior, "list")
      if (length(prior) == 1) {
        prior <- rep(prior, lng)
      } else if (length(prior) != lng) {
        stop(sprintf("Length of 'prior' must be '1' or %d.", lng))
      }
    }
  }

  params <- Map(pmcmc_parameter, name, initial, min, max, discrete, prior)
  names(params) <- populations
  class(params) <- "pmcmc_varied_parameter"

  params
}
