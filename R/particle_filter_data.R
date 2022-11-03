##' Prepare data for use with the [`particle_filter`].  This function
##' is required to use the particle filter as helps arrange data and
##' be explicit about the off-by-one errors that can occur.  It takes
##' as input your data to compare against a model, including some
##' measure of "time".  We need to convert this time into model time
##' steps (see Details).
##'
##' We require that the time variable increments in unit steps; this
##' may be relaxed in future to even steps, or possibly irregular
##' steps, but for now this assumption is required.  We assume that
##' the data in the first column is recorded at the end of a period of
##' 1 time unit.  So if you have in the first column `t = 10, data =
##' 100` we assume that the model steps from `t = 9` to to `t = 10`
##' and at that period the data has value 100.
##'
##' For continuous time models, time is simple to think about; time is
##' continuous (and real-valued) and really any time is
##' acceptable. For discrete time models there are two correlated
##' measures of time we need to consider - (1) the `dust` "time step",
##' a non-negative integer value that increases in unit steps, and (2)
##' the "model time" which is related to the dust time step based on
##' the `rate` parameter here as `<model time> = <dust time> *
##' <rate>`.  For a concrete example, consider a model where we want
##' to think in terms of days, but which we take 10 steps per
##' day. Time step 0 and model time 0 are the same, but day 1 occurs
##' at step 10, day 15 at step 150 and so on.
##'
##' @title Prepare data for use with particle filter
##'
##' @param data A [data.frame()] of data
##'
##' @param time The name of a column within `data` that represents
##'   your measure of time.  This column must be integer-like. To
##'   avoid confusion, this cannot be called `step`, `time`, or
##'   `model_time`.
##'
##' @param rate The number of model "time steps" that occur between
##'   each time point (in model time `time`).  This must also be
##'   integer-like for discrete time models and must be `NULL` for
##'   continuous time models.
##'
##' @param initial_time Optionally, an initial time to start the model
##'   from.  Provide this if you need to burn the model in, or if
##'   there is a long period with no data at the beginning of the
##'   simulation.  If provided, it must be a non-negative integer and
##'   must be at most equal to the first value of the `time` column,
##'   minus 1 (i.e., `data[[time]] - 1`). For discrete time models,
##'   this is expressed in model time.
##'
##' @param population Optionally, the name of a column within `data` that
##'   represents different populations. Must be a factor.
##'
##' @return If `population` is NULL, a data.frame with new columns
##'   `time_start` and `time_end` (required by [`particle_filter`]),
##'   along side all previous data except for the time variable, which
##'   is replaced by new `<time>_start` and `<time>_end` columns. If
##'   `population` is not NULL then a named list of data.frames as
##'   described above where each element represents populations in the
##'   order specified in the data.
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
                                 population = NULL) {
  assert_is(data, "data.frame")

  assert_scalar_character(time)
  if (!(time %in% names(data))) {
    stop(sprintf("Did not find column '%s', representing time, in data",
                 time))
  }
  if (time %in% c("step", "time", "model_time")) {
    stop(sprintf("The time column cannot be called '%s'", time))
  }

  is_continuous <- is.null(rate)
  if (is_continuous) {
    time_type <- "particle_filter_data_continuous"
  } else {
    rate <- assert_scalar_positive_integer(rate)
    time_type <- "particle_filter_data_discrete"
  }

  if (is.null(population)) {
    nesting_type <- "particle_filter_data_single"
    populations <- NULL
    model_time_end <- data[[time]]
  } else {
    nesting_type <- "particle_filter_data_nested"
    assert_scalar_character(population)
    if (!(population %in% names(data))) {
      stop(sprintf(
        "Did not find column '%s', representing population, in data",
        population))
     }
    if (!is.factor(data[[population]])) {
      stop(sprintf("Column '%s' must be a factor", population))
    }
    populations <- levels(data[[population]])

    data <- data[order(data[[population]], data[[time]]), ]
    time_split <- split(data[[time]], data[[population]])
    if (length(unique(time_split)) != 1) {
      stop("Unequal time between populations")
    }
    model_time_end <- time_split[[1]]
  }

  assert_integer(model_time_end, name = sprintf("data$%s", time))
  if (!is_continuous && !all(diff(model_time_end) == 1)) {
    ## It's possible that we can make this work ok for irregular time
    ## units, but we make this assumption below when working out the
    ## start and end step (i.e., that we assume that the data
    stop("Expected each time difference to be one unit")
  }

  if (length(model_time_end) < 2) {
    stop("Expected at least two time windows")
  }

  ## NOTE: test is against 1 because we'll start at 1 - 1 = 0
  if (!is_continuous && model_time_end[[1L]] < 1) {
    stop(sprintf("The first time must be at least 1 (but was given %d)",
                 model_time_end[[1L]]))
  }

  if (is.null(initial_time)) {
    if (is_continuous) {
      stop("'initial_time' must be given for continuous models")
    }
    initial_time <- model_time_end[[1L]] - 1
  } else {
    initial_time <- assert_integer(initial_time)
    if (initial_time < 0) {
      stop("'initial_time' must be non-negative")
    }
    if (initial_time > model_time_end[[1L]] - 1) {
      stop(sprintf("'initial_time' must be <= %d", model_time_end[[1L]] - 1))
    }
  }

  model_time_start <- c(initial_time, model_time_end[-length(model_time_end)])
  model_times <- cbind(model_time_start, model_time_end, deparse.level = 0)

  times <- model_times * (if (is_continuous) 1 else rate)
  ret <- data.frame(model_time_start = model_times[, 1],
                    model_time_end = model_times[, 2],
                    time_start = times[, 1],
                    time_end = times[, 2],
                    data[names(data) != time],
                    stringsAsFactors = FALSE, check.names = FALSE)

  names(ret)[1:2] <- paste0(time, c("_start", "_end"))
  attr(ret, "rate") <- rate
  attr(ret, "time") <- time
  attr(ret, "times") <- times
  attr(ret, "model_times") <- model_times
  attr(ret, "population") <- population
  attr(ret, "populations") <- populations
  class(ret) <- c(nesting_type, time_type, "particle_filter_data", "data.frame")

  ret
}


particle_filter_data_split <- function(data, compiled_compare) {
  population <- attr(data, "population")

  ## Drop off lots of attributes that are just annoying after the data
  ## has been split:
  attr(data, "rate") <- NULL
  attr(data, "time") <- NULL
  attr(data, "times") <- NULL
  attr(data, "model_time") <- NULL
  attr(data, "population") <- NULL
  attr(data, "populations") <- NULL
  class(data) <- "data.frame"

  if (compiled_compare) {
    dust::dust_data(data, "time_end", population)
  } else if (is.null(population)) {
    lapply(unname(split(data, seq_len(nrow(data)))), as.list)
  } else {
    ## largely copied from dust::dust_data
    rows <- lapply(seq_len(nrow(data)), function(i) as.list(data[i, ]))
    rows_grouped <- unname(split(rows, data[[population]]))
    lapply(seq_len(nrow(data) / nlevels(data[[population]])),
           function(i) lapply(rows_grouped, "[[", i))
  }
}


##' @export
`[.particle_filter_data` <- function(x, i, j, ...) { # nolint
  ret <- NextMethod("[")
  ## It's hard to detect this based on 'i' and 'j' but we can detect
  ## the effect of subsetting fairly efficiently:
  if (!identical(x$time_start, ret$time_start)) {
    ## TODO (#180): Stricter checks to come on the subset.
    k <- seq_len(nrow(x))[i]
    if (!is.null(attr(x, "population"))) {
      k <- k[k <= nrow(attr(x, "times"))]
    }
    attr(ret, "model_times") <- attr(x, "model_times")[k, , drop = FALSE]
    attr(ret, "times") <- attr(x, "times")[k, , drop = FALSE]
  }
  ret
}
