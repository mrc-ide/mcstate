context("particle_filter_data")


test_that("particle filter data validates time", {
  d <- data.frame(t = 0:10, y = 0:10)
  expect_error(
    particle_filter_data(NULL, "t", 10),
    "'data' must be a data.frame")
  expect_error(
    particle_filter_data(d, "time", 10),
    "Did not find column 'time', representing time, in data")
  expect_error(
    particle_filter_data(d, "t", 2.3),
    "'rate' must be an integer")
  expect_error(
    particle_filter_data(d + 0.5, "t", 10),
    "'t' must be an integer")
})


test_that("particle filter data creates data", {
  d <- data.frame(day = 0:10, data = 1:11)
  res <- particle_filter_data(d, "day", 10)
  expect_setequal(
    names(res),
    c("day_start", "day_end", "step_start", "step_end", "data"))
  expect_equal(res$day_start, 0:9)
  expect_equal(res$day_end, 1:10)
  expect_equal(res$step_start, 0:9 * 10)
  expect_equal(res$step_end, 1:10 * 10)
  expect_equal(res$data, 2:11)
})
