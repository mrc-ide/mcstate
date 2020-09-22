##' Thin results of running [pmcmc()]. This function may be useful
##' before using [pmcmc_predict()], or before saving pmcmc output to
##' disk.  `pmcmc_thin` takes every `thin`'th sample, while
##' `pmcmc_sample` randomly selects a total of `n_sample` samples.
##'
##' @title Thin a pmcmc chain
##' @param object Results of running [pmcmc()]
##'
##' @param burnin Optional integer number of iterations to discard as
##'   "burn-in". If given then samples `1:burnin` will be
##'   excluded from your results. Remember that the first sample
##'   represents the starting point of the chain. It is an error if
##'   this is not a positive integer or is greater than or equal to
##'   the number of samples (i.e., there must be at least one sample
##'   remaining after discarding burnin).
##'
##' @param thin Optional integer thinning factor. If given, then every
##'   `thin`'th sample is retained (e.g., if `thin` is 10
##'   then we keep samples 1, 11, 21, ...).
##'
##' @export
pmcmc_thin <- function(object, burnin = NULL, thin = NULL) {
  assert_is(object, "mcstate_pmcmc")
  i <- rep_len(TRUE, length(object$iteration))

  if (!is.null(burnin)) {
    assert_scalar_positive_integer(burnin)
    burnin_max <- max(object$iteration)
    if (burnin >= burnin_max) {
      stop(sprintf("'burnin' must be less than %d for your results",
                   burnin_max))
    }
    i <- i & object$iteration >= burnin
  }

  if (!is.null(thin)) {
    assert_scalar_positive_integer(thin)
    offset <- min(object$iteration[i])
    i <- i & ((object$iteration - offset) %% thin == 0)
  }

  pmcmc_filter(object, i)
}


##' @param n_sample The number of samples to draw from `object` *with
##'   replacement*. This means that `n_sample` can be larger than the
##'   total number of samples taken (though it probably should not)
##'
##' @export
##' @rdname pmcmc_thin
pmcmc_sample <- function(object, n_sample, burnin = NULL) {
  object <- pmcmc_thin(object, burnin)
  i <- sample.int(nrow(object$pars), n_sample, replace = TRUE)
  pmcmc_filter(object, i)
}


pmcmc_filter <- function(object, i) {
  if (!is.null(object$chain)) {
    object$chain <- object$chain[i]
  }
  object$iteration <- object$iteration[i]
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


##' Combine multiple [pmcmc()] samples into one object
##'
##' @title Combine pmcmc samples
##'
##' @param ... Arguments representing [pmcmc()] sample, i.e.,
##'   `mcstate_pmcmc` objects. Alternatively, pass a list as the
##'   argument `samples`. Names are ignored.
##'
##' @param samples A list of `mcstate_pmcmc` objects. This is often
##'   more convenient for programming against than `...`
##'
##' @export
pmcmc_combine <- function(..., samples = list(...)) {
  if (!all(vlapply(samples, inherits, "mcstate_pmcmc"))) {
    stop("All elements of '...' must be 'mcstate_pmcmc' objects")
  }
  if (any(!vlapply(samples, function(x) is.null(x$chain)))) {
    stop("Chains have already been combined")
  }
  if (length(samples) < 2) {
    stop("At least 2 samples objects must be provided")
  }
  if (length(unique(lapply(samples, function(x) colnames(x$pars)))) != 1L) {
    stop("All parameters must have the same names")
  }
  if (length(unique(vnapply(samples, function(x) nrow(x$pars)))) != 1L) {
    stop("All chains must have the same length")
  }

  iteration <- lapply(samples, "[[", "iteration")
  if (length(unique(iteration)) != 1L) {
    stop("All chains must have the same iterations")
  }
  iteration <- unlist(iteration, FALSE, FALSE)

  chain <- rep(seq_along(samples), each = nrow(samples[[1]]$pars))

  pars <- do.call(rbind, lapply(samples, "[[", "pars"))
  probabilities <- do.call(rbind, lapply(samples, "[[", "probabilities"))

  state <- lapply(samples, "[[", "state")
  if (!all_or_none(vlapply(state, is.null))) {
    stop("If 'state' is present for any samples, it must be present for all")
  }
  if (is.null(state[[1]])) {
    state <- NULL
  } else {
    state <- do.call(cbind, lapply(samples, "[[", "state"))
  }

  trajectories <- lapply(samples, "[[", "trajectories")
  if (!all_or_none(vlapply(trajectories, is.null))) {
    stop(paste("If 'trajectories' is present for any samples, it must be",
               "present for all"))
  }
  if (is.null(trajectories[[1]])) {
    trajectories <- NULL
  } else {
    trajectories <- combine_trajectories(trajectories)
  }

  ## Use the last state for predict as that will probably have most
  ## advanced seed.
  ##
  ## We might check index, rate and step here though.
  predict <- last(samples)$predict

  mcstate_pmcmc(pars, probabilities, state, trajectories, predict,
                chain, iteration)
}


combine_trajectories <- function(trajectories) {
  base <- lapply(trajectories, function(x) x[names(x) != "state"])
  if (length(unique(base)) != 1L) {
    stop("trajectories data is inconsistent")
  }

  ## This is nasty:
  trajectories_state <- lapply(trajectories, function(x)
    aperm(x$state, c(1, 3, 2)))
  trajectories_state <- array(
    unlist(trajectories_state),
    dim(trajectories_state[[1]]) * c(1, 1, length(trajectories)))
  state <- aperm(trajectories_state, c(1, 3, 2))
  rownames(state) <- rownames(trajectories[[1]]$state)

  ret <- trajectories[[1]]
  ret$state <- state
  ret
}