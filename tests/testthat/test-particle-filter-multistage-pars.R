test_that("Parameters handling works", {
  expect_error(
    particle_filter_pars_nested(list(), 2),
    "'pars' must have length 2")
  expect_equal(
    particle_filter_pars_nested(list(list(a = 1)), 1),
    list(list(a = 1)))
  expect_equal(
    particle_filter_pars_nested(list(list(a = 1), list(a = 2)), 2),
    list(list(a = 1), list(a = 2)))

  ## To test how things combine we need two of each types

  f <- function(...) identity(..1)

  ## Basic parameters
  p1a <- list(a = 1)
  p1b <- list(a = 10)

  ## Trivial multistage parameters
  p2a <- multistage_parameters(list(a = 2), list())
  p2b <- multistage_parameters(list(a = 20), list())

  ## With a single epoch changing parameters, but not transform
  p3a <- multistage_parameters(list(a = 3),
                               list(multistage_epoch(10, p1a)))
  p3b <- multistage_parameters(list(a = 30),
                               list(multistage_epoch(10, p1b)))

  ## With a single epoch chnging both parameters and transform
  p4a <- multistage_parameters(list(a = 4),
                               list(multistage_epoch(10, p1a, f)))
  p4b <- multistage_parameters(list(a = 40),
                               list(multistage_epoch(10, p1b, f)))

  ## With a single epoch chnging transform, but not parameters
  p5a <- multistage_parameters(list(a = 5),
                               list(multistage_epoch(10, NULL, f)))
  p5b <- multistage_parameters(list(a = 60),
                               list(multistage_epoch(10, NULL, f)))

  others <- c("start", "transform_state")

  res <- particle_filter_pars_nested(list(p2a), 1)
  expect_s3_class(res, "multistage_parameters")
  expect_length(res, 1)
  expect_equal(res[[1]]$pars[[1]], p2a[[1]]$pars)
  expect_equal(res[[1]][others], p2a[[1]][others])

  res <- particle_filter_pars_nested(list(p2a, p2b), 2)
  expect_s3_class(res, "multistage_parameters")
  expect_length(res, 1) # still one stage
  expect_equal(res[[1]]$pars[[1]], p2a[[1]]$pars)
  expect_equal(res[[1]]$pars[[2]], p2b[[1]]$pars)
  expect_equal(res[[1]][others], p2a[[1]][others])

  res <- particle_filter_pars_nested(list(p3a, p3b), 2)
  expect_s3_class(res, "multistage_parameters")
  expect_length(res, 2)
  expect_equal(lapply(res, function(x) x$pars),
               Map(function(a, b) list(a$pars, b$pars), p3a, p3b))
  expect_equal(lapply(res, function(x) x[others]),
               lapply(p3a, function(x) x[others]))

  res <- particle_filter_pars_nested(list(p4a, p4b), 2)
  expect_s3_class(res, "multistage_parameters")
  expect_length(res, 2)
  expect_equal(lapply(res, function(x) x$pars),
               Map(function(a, b) list(a$pars, b$pars), p4a, p4b))
  expect_equal(lapply(res, function(x) x[others]),
               lapply(p4a, function(x) x[others]))

  res <- particle_filter_pars_nested(list(p5a, p5b), 2)
  expect_s3_class(res, "multistage_parameters")
  expect_length(res, 2)
  expect_equal(lapply(res, function(x) x$pars),
               list(list(p5a[[1]]$pars, p5b[[1]]$pars), NULL))
  expect_equal(lapply(res, function(x) x[others]),
               lapply(p5a, function(x) x[others]))

  expect_error(
    particle_filter_pars_nested(list(p1a, p2a), 2),
    "'pars' must be either all multistage or all non-multistage")
  expect_error(
    particle_filter_pars_nested(list(p2a, p3a), 2),
    "Incompatible numbers of stages in pars: found 1, 2 stages")

  p_start <- multistage_parameters(list(a = 5), list(multistage_epoch(11)))
  expect_error(
    particle_filter_pars_nested(list(p3a, p_start), 2),
    "Incompatible 'start' time at phase 2")

  expect_error(
    particle_filter_pars_nested(list(p3a, p4a), 2),
    "Incompatible 'transform_state' at phase 2")

  p_no_pars <- multistage_parameters(list(a = 5), list(multistage_epoch(10)))
  expect_error(
    particle_filter_pars_nested(list(p3a, p_no_pars), 2),
    "Incompatible 'pars' at phase 2")
})
