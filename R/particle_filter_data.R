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
##' @return A data.frame with new columns `step_start` and
##'   `step_end` (required by [`particle_filter`]),
##'   along side all previous data except for the time variable, which
##'   is replaced by new `<time>_start` and `<time>_end`
##'   columns.
##'
##' @export
##' @examples
##' d <- data.frame(day = 5:20, y = runif(16))
##' mcstate::particle_filter_data(d, "day", 4)
##'
##' # If providing an initial day, then the first epoch of simulation
##' # will be longer (see the first row)
##' mcstate::particle_filter_data(d, "day", 4, 0)
particle_filter_data <- function(data, time, rate, initial_time = NULL) {
  assert_is(data, "data.frame")
  if (!(time %in% names(data))) {
    stop(sprintf("Did not find column '%s', representing time, in data",
                 time))
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

  ret <- data.frame(time_start = time_start,
                    time_end = time_end,
                    step_start = time_start * rate,
                    step_end = time_end * rate,
                    data[names(data) != time])
  names(ret)[1:2] <- paste0(time, c("_start", "_end"))
  attr(ret, "rate") <- rate
  attr(ret, "time") <- time
  class(ret) <- c("particle_filter_data", "data.frame")

  ret
}
