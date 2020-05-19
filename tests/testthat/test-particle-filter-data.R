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
