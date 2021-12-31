pfsn_initial <- function(model, initial, pars, n_particles) {
  state <- Map(initial, model$info(), rep(n_particles, length(pars)), pars)
  if (all(vlapply(state, is.null))) {
    return()
  }

  if (any(vlapply(state, is.list))) {
    stop("Setting 'step' from initial no longer supported")
  }

  len <- lengths(state)
  if (length(unique(len)) != 1) {
    stop(sprintf("initial() produced unequal state lengths %s",
                 str_collapse(len)))
  }

  if (is.null(dim(state[[1]]))) {
    state_array <- matrix(unlist(state, FALSE, FALSE), ncol = length(pars))
  } else {
    state_array <- list_to_array(state)
  }

  state_array
}


pfsn_index <- function(model, index) {
  if (is.null(index)) {
    return(NULL)
  }
  index_data <- lapply(model$info(), index)

  nok <- !all(vlapply(index_data[-1], identical, index_data[[1]]))
  if (nok) {
    stop("index must be identical across populations")
  }

  index_data[[1]]
}


## Can split and merge below here to come up with something better
## really.  All the wrangling here should come through so that compare
## does the entire thing through from comparison to weights,
## likelihood, but not kappa
pfsn_compare <- function(state, compare, data, pars) {
  log_weights <- lapply(seq_len(nlayer(state)), function(i)
    compare(array_drop(state[, , i, drop = FALSE], 3), data[[i]], pars[[i]]))
  if (all(lengths(log_weights) == 0)) {
    return(NULL)
  }

  weights <- lapply(log_weights, scale_log_weights)
  list(average = vnapply(weights, "[[", "average"),
       weights = vapply(weights, function(x) x$weights, numeric(ncol(state))))
}


## This can be used safely for the single parameter case too I think
pfsn_early_exit <- function(log_likelihood, min_log_likelihood) {
  if (any(log_likelihood == -Inf)) {
    return(TRUE)
  }
  if (length(min_log_likelihood) == 1) {
    sum(log_likelihood) < min_log_likelihood
  } else {
    all(log_likelihood < min_log_likelihood)
  }
}
