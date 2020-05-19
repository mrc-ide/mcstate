##' Prepare data for use with the \code{\link{particle_filter}}.  This
##' function is not required to use the particle filter but it helps
##' arrange data and be explicit about the off-by-one errors that can
##' occur.  It takes as input your data to compare against a model,
##' including some measure of "time".  We need to convert this time
##' into model steps.
##'
##' @title Prepare data for use with particle filter
##'
##' @param data A \code{\link{data.frame}} of data
##'
##' @param time The name of a column within \code{data} that
##'   represents your measure of time.  This column must be
##'   integer-like
##'
##' @param rate The number of model "steps" that occur between each
##'   time point (in \code{time}).  This must also be integer-like
##'
##' @export
##' @author Rich Fitzjohn
##'
##' @return A data.frame with new columns \code{step_start} and
##'   \code{step_end} (required by \code{\link{particle_filter}}),
##'   along side all previous data except for the time variable, which
##'   is replaced by new \code{<time>_start} and \code{<time>_end}
##'   columns.
##'
##' @export
##' @examples
##' d <- data.frame(day = 5:20, y = runif(16))
##' mcstate::particle_filter_data(d, "day", 4)
particle_filter_data <- function(data, time, rate) {
  assert_is(data, "data.frame")
  if (!(time %in% names(data))) {
    stop(sprintf("Did not find column '%s', representing time, in data",
                 time))
  }

  t <- data[[time]]
  assert_integer(t)
  assert_strictly_increasing(t)
  assert_integer(rate)

  time_start <- t[-nrow(data)]
  time_end <- t[-1]

  ret <- data.frame(time_start = time_start,
                    time_end = time_end,
                    step_start = time_start * rate,
                    step_end = time_end * rate,
                    data[-1, ][names(data) != time])
  names(ret)[1:2] <- paste0(time, c("_start", "_end"))
  ret
}
