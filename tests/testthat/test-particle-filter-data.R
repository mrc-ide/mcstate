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