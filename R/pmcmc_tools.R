##' Thin results of running [pmcmc()]. This function may be useful
##' before using [pmcmc_predict()], or before saving pmcmc output to
##' disk.  `pmcmc_thin` takes every `thin`'th sample, while
##' `pmcmc_sample` randomly selects a total of `n_sample` samples.
##'
##' @title Thin a pmcmc chain
##' @param object Results of running [pmcmc()]
##'
##' @param burnin Optional integer number of iterations to discard as
##'   "burn-in". If given then samples `1:burnin` will be excluded
##'   from your results. It is an error if this is not a positive
##'   integer or is greater than or equal to the number of samples
##'   (i.e., there must be at least one sample remaining after
##'   discarding burnin).
##'
##' @param thin Optional integer thinning factor. If given, then every
##'   `thin`'th sample is retained (e.g., if `thin` is 10 then we keep
##'   samples 1, 11, 21, ...).  Note that this can produce surprising
##'   results as it will always select the first sample but not
##'   necessarily always the last.
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
    i <- i & object$iteration > burnin
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
  if (is_3d_array(object$pars)) {
    i <- sample.int(nlayer(object$pars), n_sample, replace = TRUE)
  } else {
    i <- sample.int(nrow(object$pars), n_sample, replace = TRUE)
  }
  pmcmc_filter(object, i)
}


pmcmc_filter <- function(object, i) {
  if (!is.null(object$chain)) {
    object$chain <- object$chain[i]
  }
  object$iteration <- object$iteration[i]

  object$pars <- array_first_dimension(object$pars, i)
  object$probabilities <- array_first_dimension(object$probabilities, i)
  if (!is.null(object$state)) {
    object$state <- array_last_dimension(object$state, i)
  }
  if (!is.null(object$trajectories)) {
    k <- length(dim(object$trajectories$state)) - 1
    object$trajectories$state <-
      array_nth_dimension(object$trajectories$state, k, i)
  }
  if (!is.null(object$restart)) {
    k <- length(dim(object$restart$state)) - 1
    object$restart$state <- array_nth_dimension(object$restart$state, k, i)
  }

  ## This must be removed (if it was present before)
  object["pars_index"] <- list(NULL)

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
  assert_list_of(samples, "mcstate_pmcmc")

  pars <- lapply(samples, "[[", "pars_full")
  probabilities <- lapply(samples, "[[", "probabilities_full")
  iteration <- lapply(samples, "[[", "iteration")
  state <- lapply(samples, "[[", "state")
  trajectories <- lapply(samples, "[[", "trajectories")
  restart <- lapply(samples, "[[", "restart")
  adaptive <- lapply(samples, "[[", "adaptive")

  check_combine(samples, iteration, state, trajectories, restart)

  iteration <- unlist(iteration, FALSE, FALSE)

  chain <- rep(seq_along(samples), each = nrow(samples[[1]]$pars))
  pars <- array_bind(arrays = pars, dimension = 1)
  probabilities <- array_bind(arrays = probabilities, dimension = 1)

  if (is.null(state[[1]])) {
    state <- NULL
  } else {
    state <- array_bind(arrays = state)
  }

  if (is.null(trajectories[[1]])) {
    trajectories <- NULL
  } else {
    trajectories <- combine_state(trajectories)
  }

  if (is.null(restart[[1]])) {
    restart <- NULL
  } else {
    restart <- combine_state(restart)
  }

  if (is.null(adaptive[[1]])) {
    adaptive <- NULL
  } else {
    nested <- samples[[1]]$nested
    adaptive <- combine_adaptive(adaptive, nested)
  }

  ## Use the last state for predict as that will probably have most
  ## advanced seed.
  ##
  ## We might check index, rate and step here though.
  predict <- last(samples)$predict

  mcstate_pmcmc(iteration, pars, probabilities, state, trajectories,
                restart, predict, adaptive, chain)
}

check_combine <- function(samples, iteration, state, trajectories, restart) {
  if (any(!vlapply(samples, function(x) is.null(x$chain)))) {
    stop("Chains have already been combined")
  }
  if (length(samples) < 2) {
    stop("At least 2 samples objects must be provided")
  }
  if (length(unique(lapply(samples, function(x) dimnames(x$pars)))) != 1L) {
    stop("All parameters must have the same dimension names")
  }

  nok <- length(unique(viapply(samples, function(x) NLAYER(x$pars)))) != 1L ||
          length(unique(viapply(samples, function(x) nrow(x$pars)))) != 1L
  if (nok) {
    stop("All chains must have the same length")
  }
  if (length(unique(iteration)) != 1L) {
    stop("All chains must have the same iterations")
  }
  if (!all_or_none(vlapply(state, is.null))) {
    stop("If 'state' is present for any samples, it must be present for all")
  }
  if (!all_or_none(vlapply(trajectories, is.null))) {
    stop(paste(
      "If 'trajectories' is present for any samples, it must be",
      "present for all"
    ))
  }
  if (!all_or_none(vlapply(restart, is.null))) {
    stop("If 'restart' is present for any samples, it must be present for all")
  }
}


combine_state <- function(x) {
  base <- lapply(x, function(el) el[names(el) != "state"])
  if (length(unique(base)) != 1L) {
    stop(sprintf("%s data is inconsistent", deparse(substitute(x))))
  }

  state <- lapply(x, "[[", "state")
  state <- array_bind(arrays = state, dimension = length(dim(state[[1]])) - 1)

  ret <- x[[1]]
  ret$state <- state
  ret
}

combine_adaptive <- function(x, nested) {
  autocorrelation <- lapply(x, "[[", "autocorrelation")
  mean <- lapply(x, "[[", "mean")
  scaling <- lapply(x, "[[", "scaling")
  weight <- lapply(x, "[[", "weight")
  
  combine_autocorr <- function(y) {
    y_names <- rownames(y[[1]])
    ret <- array_from_list(y)
    rownames(ret) <- colnames(ret) <- y_names
    ret
  }
  
  combine_mean <- function(y) {
    y_names <- names(y[[1]])
    ret <- array_from_list(y)
    rownames(ret) <- y_names
    ret
  }
  
  if (nested) {
    populations <- names(mean[[1]]$varied)
    
    autocorr_fixed <- combine_autocorr(lapply(autocorrelation, "[[", "fixed"))
    autocorr_varied <- lapply(autocorrelation, "[[", "varied")
    autocorr_varied <-
      lapply(populations,
             function (nm) combine_autocorr(lapply(autocorr_varied, "[[", nm)))
    names(autocorr_varied) <- populations
    autocorrelation <- list(fixed = autocorr_fixed,
                            varied = autocorr_varied)
      
    mean_fixed <- combine_mean(lapply(mean, "[[", "fixed"))
    mean_varied <- lapply(mean, "[[", "varied")
    mean_varied <-
      lapply(populations,
             function (nm) combine_mean(lapply(mean_varied, "[[", nm)))
    names(mean_varied) <- populations
    mean <- list(fixed = mean_fixed,
                 varied = mean_varied)
    
    scaling_fixed <- unlist(lapply(scaling, "[[", "fixed"))
    scaling_varied <- lapply(scaling, "[[", "varied")
    scaling_varied <-
      lapply(populations,
             function (nm) unlist(lapply(scaling_varied, "[[", nm)))
    names(scaling_varied) <- populations
    scaling <- list(fixed = scaling_fixed,
                    varied = scaling_varied)
    
    
    weight_fixed <- unlist(lapply(weight, "[[", "fixed"))
    weight_varied <- lapply(weight, "[[", "varied")
    weight_varied <-
      lapply(populations,
             function (nm) unlist(lapply(weight_varied, "[[", nm)))
    names(weight_varied) <- populations
    weight <- list(fixed = weight_fixed,
                   varied = weight_varied)
    
  } else {
    autocorrelation <- combine_autocorrelation(autocorrelation)
    mean <- combine_mean(mean)
    scaling <- unlist(scaling)
    weight <- unlist(weight)
  }
  
  list(autocorrelation = autocorrelation,
       mean = mean,
       scaling = scaling,
       weight = weight)
}
