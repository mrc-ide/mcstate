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
##'   arguments: (1) the current model state, (2) the result of the
##'   `$info()` method from the model used to this point, (3) the
##'   result of the `$info()` method for the new model that was
##'   created with `pars` which will be run from this point.  Future
##'   versions of this interface may allow passing the parameters in
##'   too.
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
  ## TODO (#171): having to do this is pretty ugly, but we need both
  ## step and time below here.  Probably a little earlier processing
  ## would remove the need to do this.
  if (inherits(data, "particle_filter_data_nested")) {
    population <- data[[attr(data, "population")]]
    data <- data[population == population[[1]], ]
  }
  ## There's an awkward bit of bookkeeping here; we need to find out
  ## when each phase *ends*, as that is the index that we do the
  ## switch at.

  times <- attr(data, "model_times")

  time_pars <- vnapply(pars[-1], "[[", "start")
  time_pars_start <- c(-Inf, time_pars)
  time_pars_end <- c(time_pars, Inf)

  drop <- time_pars_end <= times[1, 1] | time_pars_start > last(times[, 2])

  if (any(drop)) {
    pars <- pars[!drop]
    time_pars <- vnapply(pars[-1], "[[", "start")
  }
  index <- which(!drop)

  time_index <- c(match(time_pars, times[, 2]), nrow(times))

  if (any(is.na(time_index))) {
    stop(sprintf("Could not map epoch to filter time: error for stage %s",
                 paste(index[is.na(time_index)], collapse = ", ")))
  }

  ## this is a bookkeeping and interpretation nightmare so disallow it:
  err <- intersect(save_restart, time_pars)
  if (length(err) > 0) {
    stop(sprintf("save_restart cannot include epoch change: error for %s",
                 paste(err, collapse = ", ")))
  }

  ## With that ruled out, the bookkeeping to split the restart dates
  ## over epochs is tolerable:
  save_restart_stage <- findInterval(save_restart, c(0, time_pars))
  save_restart_index <- seq_along(save_restart)

  for (i in seq_along(pars)) {
    pars[[i]]$time_index <- time_index[[i]]
    if (!is.null(save_restart)) {
      pars[[i]]$restart_index <- save_restart_index[save_restart_stage == i]
    }
  }

  pars
}


join_histories <- function(history, stages) {
  time_index <- vnapply(stages, "[[", "time_index")

  if (is.null(names(history[[1]]$index))) {
    ## it is possible but unlikely that we could use non-named things
    ## sensibly. We don't need that in sircovid yet so ignoring for
    ## now.
    stop("Named index required")
  }

  nms <- unique(unlist(lapply(history, function(x) names(x$index))))

  ## TODO (#172): The only use of history$index later is for names. We
  ## should refactor that to make it history_names

  ## This is not present for deterministic models!
  has_order <- !is.null(history[[1]]$order)

  value <- array(NA_real_, c(length(nms), dim(history[[1]]$value)[-1]))
  if (has_order) {
    order <- array(NA_integer_, dim(history[[1]]$order))
  }
  index <- rep(NA_integer_, length(nms))
  names(index) <- nms

  end <- time_index + 1L
  start <- c(1L, end[-length(end)] + 1L)

  rank <- length(dim(value))

  for (i in seq_along(history)) {
    j <- match(names(history[[i]]$index), nms)
    k <- seq(start[[i]], end[[i]])
    value_i <- array_last_dimension(history[[i]]$value, k)
    if (rank == 3) {
      value[j, , k] <- value_i
    } else {
      value[j, , , k] <- value_i
    }
    if (has_order) {
      order_i <- array_last_dimension(history[[i]]$order, k)
      array_last_dimension(order, k) <- order_i
    }
  }

  if (has_order) {
    list(value = value, order = order, index = index)
  } else {
    rownames(value) <- names(index)
    list(value = value, index = index)
  }
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
    array_last_dimension(restart[[i]], stages[[i]]$restart_index))
  state <- state[lengths(state) > 0]

  n <- viapply(state, nrow)
  if (length(unique(n)) > 1) {
    err <- paste(sprintf("  - %d: %s rows", seq_along(n), n), collapse = "\n")
    stop("Restart state varies in size over the simulation:\n", err)
  }

  array_bind(arrays = state)
}
