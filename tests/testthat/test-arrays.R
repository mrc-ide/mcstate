context("helpers")

test_that("array_bind", {
  m2 <- random_array(c(5, 10))
  array_bind(m2[, 1:4], m2[, 5:10])

  expect_identical(array_bind(m2[, 1:4], m2[, 5:10]), m2)

  m3 <- random_array(c(5, 7, 10))
  expect_identical(array_bind(m3[, , 1:4], m3[, , 5:10]), m3)

  m4 <- random_array(c(5, 7, 3, 10))
  expect_identical(array_bind(m4[, , , 1:4], m4[, , , 5:10]), m4)
})


test_that("preserve dimension names on merge", {
  drop_last_names <- function(m) {
    dn <- dimnames(m)
    dn[length(dn)] <- list(NULL)
    dimnames(m) <- dn
    m
  }

  m2 <- random_array(c(5, 10), TRUE)
  expect_identical(array_bind(m2[, 1:4], m2[, 5:10]), m2)
  expect_identical(array_bind(drop_last_names(m2[, 1:4]), m2[, 5:10]),
                   drop_last_names(m2))
  expect_identical(array_bind(m2[, 1:4], drop_last_names(m2[, 5:10])),
                   drop_last_names(m2))

  m3 <- random_array(c(5, 7, 10), TRUE)
  expect_identical(array_bind(m3[, , 1:4], m3[, , 5:10]), m3)

  m4 <- random_array(c(5, 7, 3, 10), TRUE)
  expect_identical(array_bind(m4[, , , 1:4], m4[, , , 5:10]), m4)
})


test_that("Can't merge incompatible arrays", {
  m4 <- random_array(c(5, 7, 3, 10))
  expect_error(
    array_bind(m4[, , , 1:4], m4[, , -1, 5:10]),
    "array 2 (dimension 3)",
    fixed = TRUE)
  expect_error(
    array_bind(m4[, , , 1:4], m4[-1, , , 1:4], m4[-1, -1, -1, 5:10]),
    "array 2 (dimension 1), array 3 (dimension 1, 2, 3)",
    fixed = TRUE)
  expect_error(
    array_bind(m4[, , , 1:4], random_array(c(5, 7, 10))),
    "Incompatible rank arrays (expected 4): array 2 (rank 3)",
    fixed = TRUE)
})


test_that("trivial case", {
  m4 <- random_array(c(5, 7, 3, 10))
  expect_identical(array_bind(m4), m4)
  expect_error(array_bind(), "Must provide at least one array")
})


test_that("reshape", {
  m4 <- random_array(c(5, 7, 10, 13))
  res <- array_reshape(m4, 3L, c(2, 5))
  expect_equal(dim(res), c(5, 7, 2, 5, 13))
  expect_equal(c(res), c(m4))

  expect_equal(res[2, 2, , ,  2], matrix(m4[2, 2, , 2], 2, 5))

  expect_error(
    array_reshape(m4, 10, c(2, 5)),
    "array only has 4 dimensions, can't update dimension 10")
  expect_error(
    array_reshape(m4, 3, c(2, 6)),
    "New dimensions (2, 6) imply dimension 3 has length 12 but found 10",
    fixed = TRUE)
})


test_that("reshape preserves dimnames where it can", {
  m4 <- random_array(c(5, 7, 10, 13), TRUE)
  res <- array_reshape(m4, 3L, c(2, 5))
  expect_equal(dim(res), c(5, 7, 2, 5, 13))
  expect_equal(dimnames(res),
               c(dimnames(m4)[1:2],
                 list(NULL, NULL),
                 dimnames(m4)[4]))

  dimnames(m4) <- NULL
  rownames(m4) <- letters[1:5]
  res <- array_reshape(m4, 3L, c(2, 5))
  expect_equal(dim(res), c(5, 7, 2, 5, 13))
  expect_equal(dimnames(res),
               list(letters[1:5], NULL, NULL, NULL, NULL))
})


test_that("drop spare dimensions", {
  m4 <- random_array(c(5, 1, 10, 1))
  res1 <- array_drop(m4, 2)
  expect_equal(c(m4), c(res1))
  expect_equal(dim(res1), c(5, 10, 1))

  res2 <- array_drop(m4, c(2, 4))
  expect_equal(c(m4), c(res2))
  expect_equal(dim(res2), c(5, 10))
})


test_that("preserve names when dropping dimensions", {
  m4 <- random_array(c(5, 1, 10, 1), TRUE)
  expect_equal(dimnames(array_drop(m4, 2)), dimnames(m4)[-2])
  expect_equal(dimnames(array_drop(m4, c(2, 4))), dimnames(m4)[-c(2, 4)])
})


test_that("Prevent impossible drops", {
  m4 <- random_array(c(5, 1, 10, 1))
  expect_error(array_drop(m4, c(2, 5)),
               "array only has 4 dimensions, can't update dimension 5")

  expect_error(array_drop(m4, 1),
               "Can't drop dimension 1 as it is length 5, not 1")
  expect_error(
    array_drop(m4, c(1, 3)),
    "Can't drop dimensions (1, 3) as they are length (5, 10), not 1",
    fixed = TRUE)
  expect_error(
    array_drop(m4, c(1, 2, 3)),
    "Can't drop dimensions (1, 3) as they are length (5, 10), not 1",
    fixed = TRUE)
})


test_that("flatten dimensions", {
  x <- array(1:24, c(2, 3, 4))
  expect_equal(array_flatten(x, 2:3),
               matrix(1:24, c(2, 12)))
  expect_equal(array_flatten(x, 1:2),
               matrix(1:24, c(6, 4)))
  expect_equal(array_flatten(x, 1:3),
               1:24)
})


test_that("Prevent impossible flattening", {
  x <- array(1:24, c(2, 3, 4))
  expect_error(array_flatten(x, 1:2 + 0.4),
               "'i' must be an integer")
  expect_error(array_flatten(x, 3:4),
               "Values of 'i' must be in [1, 3]",
               fixed = TRUE)
  expect_error(array_flatten(x, 4:6),
               "Values of 'i' must be in [1, 3]",
               fixed = TRUE)
  expect_error(array_flatten(x, 2),
               "i must be vector of at least length 2")
  expect_error(array_flatten(x, c(1, 3)),
               "All values of 'i' must be consecutive integers")
  expect_error(array_flatten(x, 3:2),
               "All values of 'i' must be consecutive integers")
  expect_error(array_flatten(x, c(2, 2, 3)),
               "All values of 'i' must be consecutive integers")
})



test_that("Can assign into last dimension of array", {
  m2 <- matrix(0, 3, 5)
  array_last_dimension(m2, 2) <- 1:3
  expect_equal(m2[, 2], 1:3)
  expect_equal(array_last_dimension(m2, 2), array(1:3, c(3, 1)))
  expect_true(all(m2[, -2] == 0))

  m3 <- array(0, c(3, 5, 7))
  array_last_dimension(m3, 2) <- 1:15
  expect_equal(m3[, , 2], matrix(1:15, 3, 5))
  expect_equal(array_last_dimension(m3, 2), array(1:15, c(3, 5, 1)))
  expect_true(all(m3[, , -2] == 0))

  m4 <- array(0, c(3, 5, 7, 11))
  array_last_dimension(m4, 2) <- 1:(3 * 5 * 7)
  expect_equal(m4[, , , 2], array(1:(3 * 5 * 7), c(3, 5, 7)))
  expect_equal(array_last_dimension(m4, 2),
               array(1:(3 * 5 * 7), c(3, 5, 7, 1)))
  expect_true(all(m4[, , , -2] == 0))

  m1 <- 1:5
  expect_error(
    array_last_dimension(m1, 2) <- 2,
    "Unexpected rank")
  expect_error(
    array_last_dimension(m1, 2),
    "Unexpected rank")
})


test_that("Can fetch the first dimension of array", {
  m2 <- matrix(1:15, 3, 5)
  expect_equal(array_first_dimension(m2, 2), m2[2, , drop = FALSE])
  expect_equal(array_first_dimension(m2, 2:3), m2[2:3, , drop = FALSE])

  m3 <- array(0, c(3, 5, 7))
  expect_equal(array_first_dimension(m3, 2), m3[2, , , drop = FALSE])
  expect_equal(array_first_dimension(m3, 2:3), m3[2:3, , , drop = FALSE])

  m4 <- array(0, c(3, 5, 7, 11))
  expect_equal(array_first_dimension(m4, 2), m4[2, , , , drop = FALSE])
  expect_equal(array_first_dimension(m4, 2:3), m4[2:3, , , , drop = FALSE])

  m1 <- 1:5
  expect_error(
    array_first_dimension(m1, 2),
    "Unexpected rank")
})


test_that("Can get an arbitrary dimension of an array", {
  m2 <- matrix(1:15, 3, 5)
  expect_equal(array_nth_dimension(m2, 2, 2), m2[, 2, drop = FALSE])
  expect_equal(array_nth_dimension(m2, 2, 2:3), m2[, 2:3, drop = FALSE])

  m3 <- array(0, c(3, 5, 7))
  expect_equal(array_nth_dimension(m3, 2, 2), m3[, 2, , drop = FALSE])
  expect_equal(array_nth_dimension(m3, 2, 2:3), m3[, 2:3, , drop = FALSE])

  m4 <- array(0, c(3, 5, 7, 11))
  expect_equal(array_nth_dimension(m4, 2, 2), m4[, 2, , , drop = FALSE])
  expect_equal(array_nth_dimension(m4, 2, 2:3), m4[, 2:3, , , drop = FALSE])

  m1 <- 1:5
  expect_error(
    array_nth_dimension(m1, 2, 2),
    "Unexpected rank")
  expect_error(
    array_nth_dimension(m2, 3, 2),
    "'k' must be in [1, 2]",
    fixed = TRUE)
})
