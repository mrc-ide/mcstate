##' Prepare data for use with the [`particle_filter`].  This function
##' is required to use the particle filter as helps arrange data and
##' be explicit about the off-by-one errors that can occur.  It takes
##' as input your data to compare against a model, including some
##' measure of "time".  We need to convert this time into model steps.
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
##'   time point (in `time`).  This must also be integer-like or NULL
##'    (for continuous time models)
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
                                 population = NULL) {
  assert_is(data, "data.frame")

  assert_scalar_character(time)
  if (!(time %in% names(data))) {
    stop(sprintf("Did not find column '%s', representing time, in data",
                 time))
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
    time_end <- data[[time]]
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
    time_end <- time_split[[1]]
  }

  assert_integer(time_end, name = sprintf("data$%s", time))
  if (!is_continuous && !all(diff(time_end) == 1)) {
    ## It's possible that we can make this work ok for irregular time
    ## units, but we make this assumption below when working out the
    ## start and end step (i.e., that we assume that the data
    stop("Expected each time difference to be one unit")
  }

  if (length(time_end) < 2) {
    stop("Expected at least two time windows")
  }

  ## NOTE: test is against 1 because we'll start at 1 - 1 = 0
  if (!is_continuous && time_end[[1L]] < 1) {
    stop(sprintf("The first time must be at least 1 (but was given %d)",
                 time_end[[1L]]))
  }

  time_start <- time_end - 1L
  if (is.null(initial_time)) {
    if (is_continuous) {
      stop("'initial_time' must be given for continuous models")
    }
  } else {
    initial_time <- assert_integer(initial_time)
    if (initial_time < 0) {
      stop("'initial_time' must be non-negative")
    }
    if (initial_time > time_start[[1L]]) {
      stop(sprintf("'initial_time' must be <= %d", time_start[[1L]]))
    }
    time_start[[1L]] <- initial_time
  }

  times <- cbind(time_start, time_end, deparse.level = 0)
  if (is_continuous) {
    steps <- NULL
    ret <- data.frame(time_start = time_start,
                      time_end = time_end,
                      time_start = time_start,
                      time_end = time_end,
                      data[names(data) != time],
                      stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    steps <- times * rate
    ret <- data.frame(time_start = time_start,
                      time_end = time_end,
                      step_start = time_start * rate,
                      step_end = time_end * rate,
                      data[names(data) != time],
                      stringsAsFactors = FALSE, check.names = FALSE)
  }
  names(ret)[1:2] <- paste0(time, c("_start", "_end"))
  attr(ret, "rate") <- rate
  attr(ret, "time") <- time
  attr(ret, "times") <- times
  attr(ret, "steps") <- steps
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
  attr(data, "steps") <- NULL
  attr(data, "population") <- NULL
  attr(data, "populations") <- NULL
  class(data) <- "data.frame"

  if (compiled_compare) {
    dust::dust_data(data, "step_end", population)
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
  if (!identical(x$step_start, ret$step_start)) {
    ## TODO (#180): Stricter checks to come on the subset.
    k <- seq_len(nrow(x))[i]
    if (!is.null(attr(x, "population"))) {
      k <- k[k <= nrow(attr(x, "steps"))]
    }
    attr(ret, "steps") <- attr(x, "steps")[k, , drop = FALSE]
    attr(ret, "times") <- attr(x, "times")[k, , drop = FALSE]
  }
  ret
}
