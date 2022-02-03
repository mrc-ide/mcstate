context("utils (assert)")

test_that("assert_is", {
  expect_error(assert_is("x", "foo"), "must be a foo")
  expect_silent(assert_is(structure("x", class = "foo"), "foo"))
})


test_that("assert_function", {
  expect_error(assert_function("x"), "must be a function")
  expect_silent(assert_function(function(x) x))
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


test_that("assert_numeric", {
  expect_silent(assert_numeric(1))
  expect_error(assert_numeric("1"), "must be a numeric")
  expect_error(assert_numeric(TRUE), "must be a numeric")
})


test_that("assert_character", {
  expect_silent(assert_character("string"))
  expect_error(assert_character(1), "must be a character")
  expect_error(assert_character(TRUE), "must be a character")
})


test_that("assert_list_of", {
  expect_error(assert_list_of(list("a"), "numeric"), "Elements of")
  expect_error(assert_list_of("a"), "must be")
  expect_equal(assert_list_of(list("a"), "character"), list("a"))
})


test_that("match_value", {
  expect_equal(match_value("aaa", c("aaa", "bbb", "ccc")), "aaa")
  expect_error(
    match_value("abc", c("aaa", "bbb", "ccc")),
    "'.+' must be one of 'aaa', 'bbb', 'ccc'")
})


test_that("check dimension names", {
  n1 <- c("a", "b")
  n2 <- c("c", "d", "e")
  n3 <- c("f", "g", "h", "i")
  arr <- array(1:(2 * 3 * 4), c(2, 3, 4), list(n1, n2, n3))
  expect_silent(assert_dimnames(arr, list(n1, n2, n3)))
  expect_error(
    assert_dimnames(arr, list(n1, c("x", "y", "z"), n3)),
    "Expected names of dimension 2 of 'arr' to match 'x', 'y', 'z'")
  expect_error(
    assert_dimnames(arr, list(A = n1, B = c("x", "y", "z"), C = n3)),
    "Expected names of dimension 2 of 'arr' to match B ('x', 'y', 'z')",
    fixed = TRUE)
  expect_error(
    assert_dimnames(arr, list(n1, NULL, n3)),
    "Expected names of dimension 2 of 'arr' to be empty")
  arr2 <- arr
  dimnames(arr2) <- list(n1, NULL, n3)
  expect_equal(assert_dimnames(arr2, list(n1, n2, n3)), arr)
  dimnames(arr2) <- list(n1, NULL, n3)
  expect_equal(assert_dimnames(arr2, list(n1, NULL, n3)), arr2)
})


test_that("check vector names", {
  nms <- c("a", "b")
  x <- set_names(1:2, nms)
  expect_silent(assert_dimnames(x, list(nms)))
  expect_error(
    assert_dimnames(x, list(c("x", "y"))),
    "Expected names of 'x' to match 'x', 'y'")
  expect_error(
    assert_dimnames(x, list(X = c("x", "y"))),
    "Expected names of 'x' to match X ('x', 'y')",
    fixed = TRUE)
  expect_error(
    assert_dimnames(x, list(NULL)),
    "Expected 'x' to have no names")
  expect_equal(assert_dimnames(unname(x), list(nms)), x)
  expect_equal(assert_dimnames(unname(x), list(NULL)), unname(x))
})
