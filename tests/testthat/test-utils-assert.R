context("utils (assert)")

test_that("assert_is", {
  expect_error(assert_is("x", "foo"), "must be a foo")
  expect_silent(assert_is(structure("x", class = "foo"), "foo"))
})


test_that("assert_integer", {
  expect_error(assert_integer(pi), "'pi' must be an integer")
  expect_identical(assert_integer(1L), 1L)
  expect_identical(assert_integer(1.0), 1L)
  expect_identical(assert_integer(1 + 1e-15), 1L)
})


test_that("assert_strictly_increasing", {
  expect_silent(assert_strictly_increasing(c(0, 1, 2)))
  expect_error(assert_strictly_increasing(c(0, 0, 1)),
               "must be strictly increasing")
  expect_error(assert_strictly_increasing(c(0, -1, -2)),
               "must be strictly increasing")
})


test_that("assert_scalar", {
  value <- NULL
  expect_silent(assert_scalar(1))
  expect_error(assert_scalar(value), "'value' must be a scalar")
  expect_error(assert_scalar(1:2), "must be a scalar")
})


test_that("assert_scalar_positive_integer", {
  expect_equal(assert_scalar_positive_integer(1L), 1L)
  expect_equal(assert_scalar_positive_integer(1000000L), 1000000L)

  value <- 0L
  expect_error(assert_scalar_positive_integer(value),
               "'value' must be at least 1")
})


test_that("assert_named", {
  expect_error(assert_named(1), "must be named")
  expect_error(assert_named(setNames(1:2, c("a", "a")), TRUE),
               "must have unique names")
  expect_silent(assert_named(setNames(1:2, c("a", "a")), FALSE))
})


test_that("assert_logical", {
  expect_error(assert_logical("one"), "must be a logical")
  expect_error(assert_logical(1), "must be a logical")
})


test_that("assert_character", {
  expect_silent(assert_character("string"))
  expect_error(assert_character(1), "must be a character")
  expect_error(assert_character(TRUE), "must be a character")
})
