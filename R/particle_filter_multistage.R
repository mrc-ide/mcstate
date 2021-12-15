##' Construct parameters for a multi-stage particle filter.
##'
##' @title Multistage filter parameters
##'
##' @param pars The parameters covering the period up to the first
##'   change in `epoch`.
##'
##' @param epochs A list of `multistage_epoch` objects corresponding
##'   to a new paramter regime starting at a new time point.
##'
##' @return An object of class `multistage_parameters`, suitable to
##'   pass through to the `run` method of [`mcstate::particle_filter`]
##'
##' @export
multistage_parameters <- function(pars, epochs) {
  ok <- vlapply(epochs, inherits, "multistage_epoch")
  if (any(!ok)) {
    stop("Expected all elements of 'epochs' to be 'multistage_epoch' objects")
  }

  sets_pars <- vlapply(epochs, function(x) is.null(x$pars))
  if (any(sets_pars) && !all(sets_pars)) {
    stop("If changing parameters, all epochs must contain 'pars'")
  }

  n_stages <- length(epochs) + 1L

  ret <- vector("list", n_stages)

  for (i in seq_len(n_stages)) {
    if (i == 1L) {
      pars_i <- pars
      transform_state_i <- NULL
      start_i <- NULL
    } else {
      pars_i <- epochs[[i - 1L]]$pars
      transform_state_i <- epochs[[i - 1L]]$transform_state
      start_i <- epochs[[i - 1L]]$start
    }

    ret[[i]] <- list(pars = pars_i,
                     start = start_i,
                     transform_state = transform_state_i)
  }

  class(ret) <- "multistage_parameters"
  ret
}


##' Describe an epoch within a [mcstate::multistage_parameters] object
##'
##' @title Multistage filter epoch
##'
##' @param start The start time, in units of time in your data set.
##'   These must correspond to time points with data.  The model will
##'   complete the step to this time point, then change parameters,
##'   then continue (so `start` represents the time point we move from
##'   with these parameters)
##'
##' @param pars Optional parameter object, replacing the model
##'   parameters at this point.  If `NULL` then the model parameters
##'   are not changed, and it is assumed that you will be changing
##'   model state via `transform_state`.
##'
##' @param transform_state Optional parameter transformation function.
##'   This could be used in two cases (1) arbitrary change to the
##'   model state (e.g., a one-off movement of state within particles
##'   at a given time point that would be otherwise awkward to code
##'   directly in your model), or (2) where you have provided `pars`
##'   and these imply a different model state size.  In this case you
##'   *must* provide `transform_state` to fill in new model state,
##'   move things around, or delete model state depending on how the
##'   state has changed.  This function will be passed three
##'   arguments: (1) the current model state, (2) the model used to
##'   run up to this point, (3) the new model that was created with
##'   `pars` which will be run from this point.  It is expected that
##'   the results of `$pars()` and `$info()` from these model objects
##'   will be useful for updating the model state.
##'
##' @export
multistage_epoch <- function(start, pars = NULL, transform_state = NULL) {
  if (is.null(transform_state)) {
    transform_state <- transform_state_identity
  } else {
    assert_function(transform_state)
  }

  ret <- list(start = start,
              pars = pars,
              transform_state = transform_state)

  class(ret) <- "multistage_epoch"
  ret
}


transform_state_identity <- function(x, ...) {
  x
}


filter_check_times <- function(pars, data, save_restart) {
  ## There's an awkward bit of bookkeeping here; we need to find out
  ## when each phase *ends*, as that is the index that we do the
  ## switch at.
  time_variable <- attr(data, "time")
  time_end_data <- data[[paste0(time_variable, "_end")]]
  time_start_pars <- vnapply(pars[-1], "[[", "start")
  step_index <- c(
    match(time_start_pars, time_end_data),
    length(time_end_data))

  if (all(is.na(step_index))) {
    ## I think that this only occurs now if there is no overlap?
    stop(sprintf("Could not map epoch to filter time: error for %s",
                 paste(which(is.na(step_index)), collapse = ", ")))
  }

  ## this is a bookkeeping and interpretation nightmare so disallow it:
  err <- intersect(save_restart, time_start_pars)
  if (length(err) > 0) {
    stop(sprintf("save_restart cannot include epoch change: error for %s",
                 paste(err, collapse = ", ")))
  }

  ## This is easy enough, but we don't need it right now and the
  ## bookkeeping is tedious
  if (length(save_restart) > 0 && any(is.na(step_index))) {
    stop("partial multistage + restart not yet supported")
  }

  ## With that ruled out, the bookkeeping to split the restart dates
  ## over epochs is tolerable:
  save_restart_stage <- findInterval(save_restart, c(0, time_start_pars))
  save_restart_index <- seq_along(save_restart)

  for (i in seq_along(pars)) {
    pars[[i]]$step_index <- step_index[[i]]
    if (!is.null(save_restart)) {
      pars[[i]]$restart_index <- save_restart_index[save_restart_stage == i]
    }
  }

  if (any(is.na(step_index))) {
    pars <- pars[!is.na(step_index)]
  }

  pars
}


join_histories <- function(history, stages) {
  step_index <- vnapply(stages, "[[", "step_index")

  if (is.null(names(history[[1]]$index))) {
    ## it is possible but unlikely that we could use non-named things
    ## sensibly. We don't need that in sircovid yet so ignoring for
    ## now.
    stop("Named index required")
  }

  nms <- unique(unlist(lapply(history, function(x) names(x$index))))

  ## TODO: The only use of history_index later is for names. We should
  ## refactor that to make it history_names

  value <- array(NA_real_, c(length(nms), dim(history[[1]]$value)[-1]))
  order <- array(NA_integer_, dim(history[[1]]$order))
  index <- rep(NA_integer_, length(nms))
  names(index) <- nms

  end <- step_index + 1L
  start <- c(1L, end[-length(end)] + 1L)

  for (i in seq_along(history)) {
    j <- match(names(history[[i]]$index), nms)
    k <- seq(start[[i]], end[[i]])
    value[j, , k] <- history[[i]]$value[, , k]
    order[, k] <- history[[i]]$order[, k]
  }

  list(value = value, order = order, index = index)
}


## Currently this does not allow saving restart state over models
## where the dimensionality changes over time (or rather, we can't
## save restart information over parts with different sizes).  We
## might be able to relax that later - simplest would be to no longer
## store a multidimensional array, perhaps nicer would be to create
## some weird wrapper object that allowed multidimensional arrays that
## differ in the size of their first dimension.
join_restart_state <- function(restart, stages) {
  state <- lapply(seq_along(stages), function(i)
    restart[[i]][, , stages[[i]]$restart_index, drop = FALSE])
  state <- state[lengths(state) > 0]

  n <- viapply(state, nrow)
  if (length(unique(n)) > 1) {
    err <- paste(sprintf("  - %d: %s rows", seq_along(n), n), collapse = "\n")
    stop("Restart state varies in size over the simulation:\n", err)
  }

  array_bind(arrays = state)
}
