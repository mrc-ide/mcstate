##' Thin results of running [pmcmc()]. This function may be useful
##' before using [pmcmc_predict()], or before saving pmcmc output to
##' disk.
##'
##' @title Thin a pmcmc chain
##' @param object Results of running [pmcmc()]
##'
##' @param burnin Optional integer number of iterations to discard as
##'   "burn-in". If given then samples `1:burnin` will be excluded
##'   from your results. Remember that the first sample represents the
##'   starting point of the chain. It is an error if this is not a
##'   positive integer or is greater than the number of samples.
##'
##' @param thin Optional integer thinning factor. If given, then every
##'   `thin`'th sample is retained (e.g., if `thin` is 10 then we keep
##'   samples 1, 11, 21, ...).
##'
##' @export
pmcmc_thin <- function(object, burnin = NULL, thin = NULL) {
  expect_is(object, "mcstate_pmcmc")
  n <- nrow(object$pars)
  i <- seq_len(n)
  if (!is.null(burnin)) {
    assert_scalar_positive_integer(burnin)
    if (burnin > n) {
      stop(sprintf("'burnin' can be at most %d for your results", n))
    }
    i <- i[-seq_len(burnin)]
  }
  if (!is.null(thin)) {
    assert_scalar_positive_integer(thin)
    i <- i[seq.int(from = 0, length.out = length(i)) %% thin == 0]
  }

  pmcmc_filter(object, i)
}


pmcmc_filter <- function(object, i) {
  object$pars <- object$pars[i, , drop = FALSE]
  object$probabilities <- object$probabilities[i, , drop = FALSE]
  if (!is.null(object$state)) {
    object$state <- object$state[, i, drop = FALSE]
  }
  if (!is.null(object$trajectories)) {
    object$trajectories <- object$trajectories[, , i, drop = FALSE]
  }
  object
}
