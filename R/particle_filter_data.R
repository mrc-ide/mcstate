##' Prepare data for use with the [`particle_filter`].  This
##' function is not required to use the particle filter but it helps
##' arrange data and be explicit about the off-by-one errors that can
##' occur.  It takes as input your data to compare against a model,
##' including some measure of "time".  We need to convert this time
##' into model steps.
##'
##' We require that the time variable increments in unit steps; this
##' may be relaxed in future to even steps, or possibly irregular
##' steps, but for now this assumption is required.  We assume that
##' the data in the first column is recorded at the end of a period of
##' 1 time unit.  So if you have in the first column `time = 10,
##' data = 100` we assume that the model steps from time 9 to to time
##' 10 and at that period the data has value 100.
##'
##' @title Prepare data for use with particle filter
##'
##' @param data A [data.frame()] of data
##'
##' @param time The name of a column within `data` that
##'   represents your measure of time.  This column must be
##'   integer-like
##'
##' @param rate The number of model "steps" that occur between each
##'   time point (in `time`).  This must also be integer-like
##'
##' @param initial_time Optionally, an initial time to start the model
##'   from.  Provide this if you need to burn the model in, or if
##'   there is a long period with no data at the beginning of the
##'   simulation.  If provided, it must be a non-negative integer and
##'   must be at most equal to the first value of the `time`
##'   column, minus 1 (i.e., `data[[time]] - 1`).
##'
##' @param population Optionally, the name of a column within `data` that
##'   represents different populations. Must be a factor.
##'
##' @param allow_unequal_times If `population` is not NULL and time-points are
##'   not equal between populations, then `TRUE` will cause the code to
##'   error, otherwise FALSE (default) will add equal time-points to all
##'   populations and fill missing data with NAs.
##'
##' @return If `population` is NULL, a data.frame with new columns
##'   `step_start` and `step_end` (required by [`particle_filter`]),
##'   along side all previous data except for the time variable, which
##'   is replaced by new `<time>_start` and `<time>_end`
##'   columns. If `population` is not NULL then a named list of data.frames as
##'   described above where each element represents populations in the order
##'   specified in the data.
##'
##' @export
##' @examples
##' d <- data.frame(day = 5:20, y = runif(16))
##' mcstate::particle_filter_data(d, "day", 4)
##'
##' # If providing an initial day, then the first epoch of simulation
##' # will be longer (see the first row)
##' mcstate::particle_filter_data(d, "day", 4, 0)
##'
##' # If including populations:
##' d <- data.frame(day = 5:20, y = runif(16),
##'                 population = factor(rep(letters[1:2], each = 16)))
##' mcstate::particle_filter_data(d, "day", 4, 0, "population")
particle_filter_data <- function(data, time, rate, initial_time = NULL,
                                 population = NULL,
                                 allow_unequal_times = FALSE) {

  assert_is(data, "data.frame")
  if (!(time %in% names(data))) {
    stop(sprintf("Did not find column '%s', representing time, in data",
                 time))
  }

  if (!is.null(population)) {
    if (!(population %in% names(data))) {
      stop(sprintf("Did not find column '%s', representing population, in
                   data", population))
    }
    return(particle_filter_data_nested(data, time, rate, initial_time,
                                       population, allow_unequal_times))
  }


  t <- assert_integer(data[[time]])
  if (!all(diff(t) == 1)) {
    ## It's possible that we can make this work ok for irregular time
    ## units, but we make this assumption below when working out the
    ## start and end step (i.e., that we assume that the data
    stop("Expected each time difference to be one unit")
  }
  if (t[[1L]] < 1) {
    stop(sprintf(
      "The first time must be at least 1 (but was given %d)", t[[1L]]))
  }

  if (nrow(data) < 2) {
    stop("Expected at least two time windows")
  }

  rate <- assert_scalar_positive_integer(rate)
  t <- clean_pf_times(t, initial_time, data)

  ret <- data.frame(time_start = t$time_start,
                    time_end = t$time_end,
                    step_start = t$time_start * rate,
                    step_end = t$time_end * rate,
                    data[names(data) != time])
  names(ret)[1:2] <- paste0(time, c("_start", "_end"))
  attr(ret, "rate") <- rate
  attr(ret, "time") <- time
  class(ret) <- c("particle_filter_data_single", "particle_filter_data",
                  "data.frame")

  ret
}

## TODO (#171): this needs simplification and general tidying, all a
## bit complicated still throughout here.
particle_filter_data_nested <- function(data, time, rate, initial_time,
                                        population, allow_unequal_times) {
  # catch invalid user call
  if (is.null(population)) {
    stop("'population' must be non-NULL")
  }

  if (!is.factor(data[[population]])) {
    stop(sprintf("Column '%s' must be a factor", population))
  }

  groups <- unique(data[[population]])

  assert_logical(allow_unequal_times)

  split_data <- split(data[names(data) != population], data[[population]])
  vars <- names(data)[!(names(data) %in% c(time, population))]

  nr <- viapply(split_data, nrow)
  if (length(unique(nr)) > 1) {
    if (allow_unequal_times) {
      ## assumed quicker handling one matrix then re-splitting than lapply
      times <- as.integer(sort(unique(unlist(lapply(split_data, "[[", time)))))
      split_data <- lapply(split_data, function(x) {
        df <- data.frame(matrix(ncol = ncol(x) - 1, nrow = length(times)))
        colnames(df) <- colnames(x)[-ncol(x)]
        df[[time]] <- times
        df[match(x[[time]], times), vars] <- as.matrix(x[, vars])
        df
      })
    } else {
      stop("Unequal time between populations")
    }
  }

  rate <- assert_scalar_positive_integer(rate)
  ## times now equal so select first arbitrarily
  t <- clean_pf_times(split_data[[1]][[time]], initial_time, data)

  data <- do.call(rbind, split_data)

  ## TODO (#171): this is not ideal because 'population' is then
  ## silently unavailable (a bit like step is, but worse).
  ret <- data.frame(time_start = t$time_start,
                    time_end = t$time_end,
                    step_start = t$time_start * rate,
                    step_end = t$time_end * rate,
                    population = rep(groups, each = length(t$time_start)),
                    data[vars],
                    stringsAsFactors = FALSE)
  rownames(ret) <- NULL
  names(ret)[1:2] <- paste0(time, c("_start", "_end"))
  attr(ret, "rate") <- rate
  attr(ret, "time") <- time
  ## TODO (#171): these will become flexible soon, but this at least
  ## establishes the contract for things that use the output of this
  ## function.
  attr(ret, "population") <- "population"
  attr(ret, "n_populations") <- nlevels(ret$population)
  class(ret) <- c("particle_filter_data_nested", "particle_filter_data",
                  "data.frame")
  ret
}

clean_pf_times <- function(times, initial_time, data) {
  t <- assert_integer(times)

  if (!all(diff(t) == 1)) {
    ## It's possible that we can make this work ok for irregular time
    ## units, but we make this assumption below when working out the
    ## start and end step (i.e., that we assume that the data
    stop("Expected each time difference to be one unit")
  }
  if (t[[1L]] < 1) {
    stop(sprintf(
      "The first time must be at least 1 (but was given %d)", t[[1L]]))
  }

  if (nrow(data) < 2) {
    stop("Expected at least two time windows")
  }

  time_end <- t
  time_start <- t - 1L
  if (!is.null(initial_time)) {
    initial_time <- assert_integer(initial_time)
    if (initial_time < 0) {
      stop("'initial_time' must be non-negative")
    }
    if (initial_time > time_start[[1L]]) {
      stop(sprintf(
        "'initial_time' must be less than %d", time_start[[1L]]))
    }
    time_start[[1L]] <- initial_time
  }

  list(time_start = time_start, time_end = time_end)
}
