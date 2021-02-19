context("particle_filter_data")


test_that("particle filter data validates time", {
  d <- data.frame(t = 1:11, y = 0:10)
  expect_error(
    particle_filter_data(NULL, "t", 10),
    "'data' must be a data.frame")
  expect_error(
    particle_filter_data(d, "time", 10),
    "Did not find column 'time', representing time, in data")
  expect_error(
    particle_filter_data(d + 0.5, "t", 10),
    "'data[[time]]' must be an integer",
    fixed = TRUE)
  expect_error(
    particle_filter_data(d - 1, "t", 10),
    "The first time must be at least 1 (but was given 0)",
    fixed = TRUE)
  expect_error(
    particle_filter_data(d * 2, "t", 10),
    "Expected each time difference to be one unit")
})


test_that("particle filter data validates rate", {
  d <- data.frame(t = 1:11, y = 0:10)
  expect_error(
    particle_filter_data(d, "t", 2.3),
    "'rate' must be an integer")
})


test_that("particle filter data validates initial_time", {
  d <- data.frame(t = 1:11, y = 0:10)
  expect_error(
    particle_filter_data(d, "t", 2, -10),
    "'initial_time' must be non-negative")
  expect_error(
    particle_filter_data(d, "t", 2, 2),
    "'initial_time' must be less than 0")
  expect_error(
    particle_filter_data(d, "t", 2, 0.5),
    "'initial_time' must be an integer")
})


test_that("particle filter data creates data", {
  d <- data.frame(day = 1:11, data = seq(0, 1, by = 0.1))
  res <- particle_filter_data(d, "day", 10)
  expect_setequal(
    names(res),
    c("day_start", "day_end", "step_start", "step_end", "data"))
  expect_equal(res$day_start, 0:10)
  expect_equal(res$day_end, 1:11)
  expect_equal(res$step_start, 0:10 * 10)
  expect_equal(res$step_end, 1:11 * 10)
  expect_equal(res$data, d$data)
})


test_that("particle filter can offset initial data", {
  d <- data.frame(hour = 11:20, a = runif(10), b = runif(10))
  res <- particle_filter_data(d, "hour", 4, 1)
  cmp <- data.frame(hour_start = c(1, 11:19),
                    hour_end = 11:20,
                    step_start = c(4, 11:19 * 4),
                    step_end = 11:20 * 4,
                    a = d$a,
                    b = d$b)
  attr(cmp, "rate") <- 4
  attr(cmp, "time") <- "hour"
  class(cmp) <- c("particle_filter_data", "data.frame")
  expect_equal(res, cmp)
})


test_that("require more than one observation", {
  d <- data.frame(hour = 1:2, a = 2:3, b = 3:4)
  expect_error(
    particle_filter_data(d[1, ], "hour", 10),
    "Expected at least two time windows")
  expect_silent(
    particle_filter_data(d, "hour", 10))
})

test_that("particle filter data with populations creates data - equal", {
  d <- data.frame(day = 1:11, data = seq(0, 1, by = 0.1),
                  population = rep(letters[1:2], each = 11),
                  stringsAsFactors = TRUE)
  res <- particle_filter_data(d, "day", 10, population = "population")

  expect_equal(names(res), letters[1:2])
  expect_s3_class(res, "particle_filter_data_nested")
  expect_s3_class(res[[1]], "particle_filter_data")
  expect_s3_class(res[[2]], "particle_filter_data")

  expect_setequal(
    names(res[[1]]),
    c("day_start", "day_end", "step_start", "step_end", "data"))
  expect_equal(res[[1]]$day_start, 0:10)
  expect_equal(res[[1]]$day_end, 1:11)
  expect_equal(res[[1]]$step_start, 0:10 * 10)
  expect_equal(res[[1]]$step_end, 1:11 * 10)
  expect_equal(res[[1]]$data, d$data[1:11])

  expect_setequal(
    names(res[[2]]),
    c("day_start", "day_end", "step_start", "step_end", "data"))
  expect_equal(res[[2]]$day_start, 0:10)
  expect_equal(res[[2]]$day_end, 1:11)
  expect_equal(res[[2]]$step_start, 0:10 * 10)
  expect_equal(res[[2]]$step_end, 1:11 * 10)
  expect_equal(res[[2]]$data, d$data[1:11])
})

test_that("particle filter data with populations creates data - unequal", {
  d <- data.frame(day = c(seq.int(1, 10, 2), 1:10),
                  data1 = c(runif(5), runif(10)),
                  data2 = c(runif(5), runif(10)),
                  population = rep(letters[1:2], times = c(5, 10)),
                  stringsAsFactors = TRUE)
  expect_error(particle_filter_data(d, "day", 10, population = "population",
                                    error_on_unequal = TRUE), "Unequal")

  res <- particle_filter_data(d, "day", 10, population = "population",
                              error_on_unequal = FALSE)

  expect_s3_class(res, "particle_filter_data_nested")
  expect_s3_class(res[[1]], "particle_filter_data")
  expect_s3_class(res[[2]], "particle_filter_data")

  expect_setequal(
    names(res[[1]]),
    c("day_start", "day_end", "step_start", "step_end", "data1", "data2"))
  expect_equal(res[[1]]$day_start, 0:9)
  expect_equal(res[[1]]$day_end, 1:10)
  expect_equal(res[[1]]$step_start, 0:9 * 10)
  expect_equal(res[[1]]$step_end, 1:10 * 10)
  expect_equal(subset(res[[1]], day_end %in% d$day[1:5])[, "data1"],
               d$data1[1:5])
  expect_equal(subset(res[[1]], !(day_end %in% d$day[1:5]))[, "data1"],
               rep(NA_integer_, 5))
  expect_equal(subset(res[[1]], day_end %in% d$day[1:5])[, "data2"],
               d$data2[1:5])
  expect_equal(subset(res[[1]], !(day_end %in% d$day[1:5]))[, "data2"],
               rep(NA_integer_, 5))

  expect_setequal(
    names(res[[2]]),
    c("day_start", "day_end", "step_start", "step_end", "data1", "data2"))
  expect_equal(res[[2]]$day_start, 0:9)
  expect_equal(res[[2]]$day_end, 1:10)
  expect_equal(res[[2]]$step_start, 0:9 * 10)
  expect_equal(res[[2]]$step_end, 1:10 * 10)
  expect_equal(res[[2]]$data1, d$data1[6:15])
  expect_equal(res[[2]]$data2, d$data2[6:15])
})

test_that("particle_filter_data_multi - error on NULL", {
  expect_error(particle_filter_data_multi(population = NULL), "must be non")
})