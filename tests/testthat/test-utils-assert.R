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
  x <- NULL
  expect_error(assert_scalar(x), "'x' must be a scalar")
  expect_error(assert_scalar(1:5), "must be a scalar")
  expect_silent(assert_scalar(1))
  expect_silent(assert_scalar_integer(1))
})
