##' Thin results of running \code{\link{pmcmc}}. This function may be
##' useful before using \code{\link{pmcmc_predict}}, or before saving
##' pmcmc output to disk.
##'
##' @title Thin a pmcmc chain
##' @param object Results of running \code{\link{pmcmc}}
##'
##' @param burnin Optional integer number of iterations to discard as
##'   "burn-in". If given then samples \code{1:burnin} will be
##'   excluded from your results. Remember that the first sample
##'   represents the starting point of the chain. It is an error if
##'   this is not a positive integer or is greater than or equal to
##'   the number of samples (i.e., there must be at least one sample
##'   remaining after discarding burnin).
##'
##' @param thin Optional integer thinning factor. If given, then every
##'   \code{thin}'th sample is retained (e.g., if \code{thin} is 10
##'   then we keep samples 1, 11, 21, ...).
##'
##' @export
pmcmc_thin <- function(object, burnin = NULL, thin = NULL) {
  assert_is(object, "mcstate_pmcmc")
  n <- nrow(object$pars)
  i <- seq_len(n)
  if (!is.null(burnin)) {
    assert_scalar_positive_integer(burnin)
    if (burnin >= n) {
      stop(sprintf("'burnin' can be at most %d for your results", n - 1))
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
    object$trajectories$state <- object$trajectories$state[, i, , drop = FALSE]
  }
  object
}
