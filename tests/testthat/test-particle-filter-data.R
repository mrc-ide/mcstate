context("particle_filter_data")


test_that("particle filter data validates time", {
  d <- data.frame(t = 1:11, y = 0:10)
  expect_error(
    particle_filter_data(NULL, "t", 10, 0),
    "'data' must be a data.frame")
  expect_error(
    particle_filter_data(d, "time", 10, 0),
    "Did not find column 'time', representing time, in data")
  expect_error(
    particle_filter_data(d + 0.5, "t", 10, 0),
    "'data$t' must be an integer",
    fixed = TRUE)
  expect_error(
    particle_filter_data(d - 2, "t", 10, 0),
    "All times must be non-negative",
    fixed = TRUE)
  expect_error(
    particle_filter_data(d, "t", 10, -1),
    "'initial_time' must be non-negative",
    fixed = TRUE)
  expect_error(
    particle_filter_data(d, "t", 10, 2),
    "'initial_time' must be <= 1",
    fixed = TRUE)
})


test_that("can't use reserved names for time column", {
  expect_error(
    particle_filter_data(data_frame(time = 1:10), "time", 1, 0),
    "The time column cannot be called 'time'")
  expect_error(
    particle_filter_data(data_frame(step = 1:10), "step", 1, 0),
    "The time column cannot be called 'step'")
  expect_error(
    particle_filter_data(data_frame(model_time = 1:10), "model_time", 1, 0),
    "The time column cannot be called 'model_time'")
})


test_that("particle filter data validates rate", {
  d <- data.frame(t = 1:11, y = 0:10)
  expect_error(
    particle_filter_data(d, "t", 2.3, 0),
    "'rate' must be an integer")
})


test_that("particle filter data validates initial_time", {
  d <- data.frame(t = 1:11, y = 0:10)
  expect_error(
    particle_filter_data(d, "t", 2, -10),
    "'initial_time' must be non-negative")
  expect_error(
    particle_filter_data(d, "t", 2, 2),
    "'initial_time' must be <= 1")
  expect_error(
    particle_filter_data(d, "t", 2, 0.5),
    "'initial_time' must be an integer")
})


test_that("particle filter data creates data", {
  d <- data.frame(day = 1:11, data = seq(0, 1, by = 0.1))
  res <- particle_filter_data(d, "day", 10, 0)
  expect_setequal(
    names(res),
    c("day_start", "day_end", "time_start", "time_end", "data"))
  expect_equal(res$day_start, 0:10)
  expect_equal(res$day_end, 1:11)
  expect_equal(res$time_start, 0:10 * 10)
  expect_equal(res$time_end, 1:11 * 10)
  expect_equal(res$data, d$data)
  expect_equal(attr(res, "rate"), 10)
  expect_equal(attr(res, "time"), "day")
  expect_equal(attr(res, "model_times"), cbind(0:10, 1:11, deparse.level = 0))
  expect_equal(attr(res, "times"), attr(res, "model_times") * 10)
  expect_s3_class(
    res,
    c("particle_filter_data_single",
      "particle_filter_data_discrete",
      "particle_filter_data",
      "data.frame"),
    exact = TRUE)
})


test_that("particle filter can offset initial data", {
  d <- data.frame(hour = 11:20, a = runif(10), b = runif(10))
  res <- particle_filter_data(d, "hour", 4, 1)

  cmp <- data.frame(hour_start = c(1, 11:19),
                    hour_end = 11:20,
                    time_start = c(4, 11:19 * 4),
                    time_end = 11:20 * 4,
                    a = d$a,
                    b = d$b)
  attr(cmp, "rate") <- 4
  attr(cmp, "time") <- "hour"
  attr(cmp, "model_times") <- cbind(c(1, 11:19), 11:20)
  attr(cmp, "times") <- attr(cmp, "model_times") * 4
  class(cmp) <- c("particle_filter_data_single",
                  "particle_filter_data_discrete",
                  "particle_filter_data",
                  "data.frame")
  expect_equal(res, cmp)
})


test_that("require more than one observation", {
  d <- data.frame(hour = 1:2, a = 2:3, b = 3:4)
  expect_error(
    particle_filter_data(d[1, ], "hour", 10, 0),
    "Expected at least two time windows")
  expect_silent(
    particle_filter_data(d, "hour", 10, 0))
})


test_that("particle filter data with populations creates data - equal", {
  data <- runif(22)
  d <- data.frame(day = 1:11,
                  data = data,
                  group = rep(letters[1:2], each = 11),
                  stringsAsFactors = TRUE)
  d <- d[sample.int(nrow(d)), ]
  res <- particle_filter_data(d, "day", 10, 0, population = "group")

  expect_s3_class(res, "particle_filter_data_nested")

  expect_setequal(
    names(res),
    c("day_start", "day_end", "time_start", "time_end", "group", "data"))
  expect_equal(res$day_start, rep(0:10, 2))
  expect_equal(res$day_end, rep(1:11, 2))
  expect_equal(res$time_start, rep(0:10, 2) * 10)
  expect_equal(res$time_end, rep(1:11, 2) * 10)
  expect_equal(res$group, factor(rep(c("a", "b"), each = 11)))
  expect_equal(res$data, data)

  expect_equal(attr(res, "rate"), 10)
  expect_equal(attr(res, "time"), "day")
  expect_equal(attr(res, "model_times"), cbind(0:10, 1:11))
  expect_equal(attr(res, "times"), cbind(0:10, 1:11) * 10)
  expect_equal(attr(res, "population"), "group")
  expect_equal(attr(res, "populations"), c("a", "b"))
})


test_that("particle filter data with populations creates data - unequal", {
  d <- data.frame(day = c(seq.int(1, 10, 2), 1:10),
                  data1 = c(runif(5), runif(10)),
                  data2 = c(runif(5), runif(10)),
                  population = rep(letters[1:2], times = c(5, 10)),
                  stringsAsFactors = TRUE)
  expect_error(
    particle_filter_data(d, "day", 10, population = "population"),
    "Unequal time between populations")
})


test_that("particle_filter_data_multi - errors", {
  d <- data.frame(day = 1:11,
                  data = runif(22),
                  group = rep(letters[1:2], each = 11),
                  stringsAsFactors = FALSE)
  expect_error(
    particle_filter_data(d, "day", 1, population = "grp"),
    "Did not find column 'grp', representing population, in data")
  expect_error(
    particle_filter_data(d, "day", 1, population = "group"),
    "Column 'group' must be a factor")
})


test_that("particle_filter_data for continuous time", {
  d <- data.frame(month = 4:24,
                  data = runif(21),
                  stringsAsFactors = FALSE)

  res <- particle_filter_data(d, "month", NULL, initial_time = 0)
  expect_setequal(
    names(res),
    c("month_start", "month_end", "time_start", "time_end", "data"))
  expect_equal(res$month_start, c(0, 4:23))
  expect_equal(res$month_end, 4:24)
  expect_equal(res$time_start, res$month_start)
  expect_equal(res$time_end, res$month_end)
  expect_equal(res$data, d$data)
  expect_equal(attr(res, "rate"), NULL)
  expect_equal(attr(res, "time"), "month")
  expect_equal(attr(res, "model_times"),
               cbind(c(0, 4:23), 4:24, deparse.level = 0))
  expect_s3_class(
    res,
    c("particle_filter_data_single",
      "particle_filter_data_continuous",
      "particle_filter_data",
      "data.frame"),
    exact = TRUE)
})


test_that("particle_filter_data for continuous time by month", {
  d <- data.frame(day = seq(30, by = 30, length.out = 10),
                  data = runif(10),
                  stringsAsFactors = FALSE)

  res <- particle_filter_data(d, "day", NULL, initial_time = 0)

  expect_equal(res$day_start, d$day - 30)
  expect_equal(res$day_end, d$day)
  expect_equal(res$time_start, res$day_start)
  expect_equal(res$time_end, res$day_end)
})

test_that("particle_filter_data for continuous time requires initial time", {
  d <- data.frame(month = 4:24,
                  data = runif(21),
                  stringsAsFactors = FALSE)
  expect_error(particle_filter_data(d, "month", NULL),
               "'initial_time' must be given for continuous models")
})


test_that("particle filter data can construct with non-unit time data", {
  dat <- example_sir()

  d1 <- dat$data_raw
  d1$incidence[rep(c(TRUE, FALSE), length.out = nrow(d1))] <- NA
  d2 <- d1[!is.na(d1$incidence), ]

  df1 <- particle_filter_data(d1, "day", 4, 0)
  df2 <- particle_filter_data(d2, "day", 4, 0)

  i <- which(!is.na(d1$incidence))
  expect_equal(df2$day_start, df1$day_start[i - 1])
  expect_equal(df2$day_end, df1$day_end[i])
  expect_equal(df2$time_start, df1$time_start[i - 1])
  expect_equal(df2$time_end, df1$time_end[i])
  expect_equal(df2$incidence, df1$incidence[i])

  expect_equal(attr(df2, "rate"), attr(df1, "rate"))
  expect_equal(attr(df2, "time"), attr(df1, "time"))
})

test_that("particle filter data can construct with irregular time data", {
  dat <- example_sir()

  set.seed(1)
  d1 <- dat$data_raw
  d1$incidence[c(runif(nrow(d1) - 1) < 0.5, FALSE)] <- NA
  d2 <- d1[!is.na(d1$incidence), ]

  df1 <- particle_filter_data(d1, "day", 4, 0)
  df2 <- particle_filter_data(d2, "day", 4, 0)

  i <- which(!is.na(d1$incidence))
  ## expect_equal(df2$day_start, df1$day_start[i - 1])
  expect_equal(df2$day_start, c(0, df2$day_end[-nrow(df2)]))
  expect_equal(df2$day_end, df1$day_end[i])
  expect_equal(df2$time_start, c(0, df2$time_end[-nrow(df2)]))
  expect_equal(df2$time_end, df1$time_end[i])
  expect_equal(df2$incidence, df1$incidence[i])

  expect_equal(attr(df2, "rate"), attr(df1, "rate"))
  expect_equal(attr(df2, "time"), attr(df1, "time"))
})


test_that("particle filter data warns if initial time not given", {
  d <- data.frame(day = 2:12, data = seq(0, 1, by = 0.1))
  expect_warning(
    res <- particle_filter_data(d, "day", 10),
    "'initial_time' should be provided. I'm assuming '1'")
  expect_equal(res, particle_filter_data(d, "day", 10, 1))
})
