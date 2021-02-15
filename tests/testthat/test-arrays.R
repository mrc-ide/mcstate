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
