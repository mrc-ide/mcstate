##' Bind a number of arrays, usually by their last dimension. This is
##' useful for binding together the sorts of arrays produced by dust
##' and mcstate's simulation functions.
##'
##' @title Bind arrays
##'
##' @param ... Any number of arrays. All dimensions must the the same,
##'   except for the dimension being bound on which may vary.
##'
##' @param arrays As an alternative to using `...` you can provide a
##'   list directly. This is often nicer to program with.
##'
##' @param dimension The dimension to bind on; by default `NULL` means
##'   the last dimension.
##'
##' @return A single array object
##' @export
##' @examples
##' # Consider two matricies; this is equivalent to rbind and is
##' # pretty trivial
##' m1 <- matrix(1, 4, 5)
##' m2 <- matrix(2, 4, 2)
##' mcstate::array_bind(m1, m2)
##'
##' # For a 4d array though it's less obvious
##' a1 <- array(1, c(2, 3, 4, 5))
##' a2 <- array(2, c(2, 3, 4, 1))
##' a3 <- array(3, c(2, 3, 4, 3))
##' dim(mcstate::array_bind(a1, a2, a3))
array_bind <- function(..., arrays = list(...), dimension = NULL) {
  if (length(arrays) == 0) {
    stop("Must provide at least one array")
  }
  if (length(arrays) == 1L) {
    return(arrays[[1L]])
  }

  d <- lapply(arrays, dim)
  r <- lengths(d)
  err <- r != r[[1]]
  if (any(err)) {
    detail <- sprintf("array %s (rank %d)", which(err), r[err])
    stop(sprintf("Incompatible rank arrays (expected %d): %s",
                 r[[1]], detail))
  }
  r <- r[[1]]
  if (is.null(dimension)) {
    dimension <- r
  }

  d <- matrix(unlist(d), r)

  d_shared <- d[-dimension, , drop = FALSE]
  d_diff <- d_shared - d_shared[, 1] != 0
  if (any(d_diff)) {
    err <- which(d_diff, TRUE)
    err_dim <- unname(seq_len(r)[-dimension][err[, 1]])
    err_arr <- err[, 2]
    err <- vcapply(split(err_dim, err_arr), paste, collapse = ", ")
    detail <- paste(sprintf("array %s (dimension %s)", names(err), unname(err)),
                    collapse = ", ")
    stop(paste("Incompatible dimension arrays:", detail))
  }

  n <- d[dimension, ]
  d <- d[, 1]
  d[[dimension]] <- sum(n)

  if (dimension == r) {
    ret <- array(unlist(arrays), d)
  } else {
    offset <- c(0, cumsum(n))
    ret <- array(arrays[[1]][[1]] * NA, d)
    idx <- array(seq_along(ret), d)
    for (i in seq_along(arrays)) {
      j <- seq_len(n[[i]]) + offset[[i]]
      k <- c(array_nth_dimension(idx, dimension, j))
      ret[k] <- c(arrays[[i]])
    }
  }

  dn <- dimnames(arrays[[1]])
  if (!is.null(dn)) {
    dn_this <- lapply(arrays, function(x) dimnames(x)[[dimension]])
    if (all(lengths(dn_this) > 0)) {
      dn[[dimension]] <- unlist(dn_this)
    } else {
      dn[dimension] <- list(NULL)
    }
    dimnames(ret) <- dn
  }

  ret
}


##' Reshape one dimension of a multidimensional array. Use this to say
##' that some dimension (say with length 20) actually represents a
##' number of other dimensions (e.g., 2 x 10 or 2 x 2 x 5). This might
##' be the case if you've been doing a simulation with a large number
##' of parameter sets that are pooled over some other grouping factors
##' (e.g., in a sensitivity analysis)
##'
##' @title Rehape an array dimension
##'
##' @param x An array
##'
##' @param i The index of the dimension to expand
##'
##' @param d The new dimensions for data in the i'th dimension of x
##'
##' @return A multidimensional array
##' @export
##' @seealso [mcstate::array_flatten] which undoes this operation
##' @examples
##' # Suppose we had a 4 x 6 array of data:
##' m <- matrix(1:24, 4, 6)
##'
##' # And suppose that the second dimension really represented a 2 x 3
##' # matrix; so that looking at one copy of the 2nd dimension we see
##' m[1, ]
##'
##' # But instead we might want to see
##' res <- mcstate::array_reshape(m, 2, c(2, 3))
##' res[1, , ]
array_reshape <- function(x, i, d) {
  dx <- dim(x)
  if (length(dx) < i) {
    stop(sprintf(
      "array only has %d dimensions, can't update dimension %d",
      length(dx), i))
  }
  if (dx[[i]] != prod(d)) {
    stop(sprintf(
      "New dimensions (%s) imply dimension %d has length %d but found %d",
      paste(d, collapse = ", "), i, prod(d), dx[[i]]))
  }

  dn <- dimnames(x)

  ## The actual reshape is easy:
  dim(x) <- append(dx[-i], d, i - 1L)

  ## Can't preserve dimension names on modified dimensions
  if (!is.null(dn)) {
    dimnames(x) <- append(dn[-i], rep(list(NULL), length(d)), i - 1L)
  }
  x
}


##' Drop specific array dimensions that are equal to 1. This a more
##' explicit, safer version of [drop], which requires you indicate
##' which dimensions will be dropped and errors if dimensions can't be
##' dropped.
##'
##' @title Drop specific array dimensions
##'
##' @param x An array
##'
##' @param i Index or indices of dimensions to remove
##'
##' @return An array
##' @export
##' @examples
##'
##' # Suppose we have an array with a redundant 2nd dimension
##' m <- array(1:25, c(5, 1, 5))
##'
##' # commonly we might drop this with
##' drop(m)
##'
##' # in this case, array_drop is the same:
##' mcstate::array_drop(m, 2)
##'
##' # However, suppose that our matrix had, in this case, a first
##' # dimension that was also 1 but we did not want to drop it:
##' m2 <- m[1, , , drop = FALSE]
##'
##' # Here, drop(m2) returns just a vector, discarding our first dimension
##' drop(m2)
##'
##' # However, array_drop will preserve that dimension
##' mcstate::array_drop(m2, 2)
array_drop <- function(x, i) {
  dx <- dim(x)
  rank <- length(dx)
  if (rank < max(i)) {
    stop(sprintf(
      "array only has %d dimensions, can't update dimension %d",
      rank, max(i)))
  }
  if (!all(dx[i] == 1)) {
    err <- i[dx[i] != 1]
    if (length(err) == 1L) {
      stop(sprintf(
        "Can't drop dimension %d as it is length %d, not 1",
        err, dx[err]))
    } else {
      stop(sprintf(
        "Can't drop dimensions (%s) as they are length (%s), not 1",
        paste(err, collapse = ", "), paste(dx[err], collapse = ", ")))
    }
  }
  dn <- dimnames(x)
  dim(x) <- dim(x)[-i]
  if (!is.null(dn)) {
    if (rank == 2) {
      x <- c(x)
      names(x) <- dn[[-i]]
    } else {
      dimnames(x) <- dn[-i]
    }
  }
  x
}


##' Flatten array dimensions into a single dimension. This takes a
##' multidimensional array and converts some dimensions of it into a
##' vector. Use this to drop out "middle" dimensions of a structured
##' array. This is conceptually the inverse of [mcstate::array_reshape]
##'
##' @title Flatten array dimensions
##' @param x An array
##'
##' @param i An integer vector of dimensions to flatten
##'
##' @return A new array with at one or more dimensions removed
##' @seealso [mcstate::array_flatten] which adds structure
##' @export
##' @examples
##' x <- array(1:12, c(2, 3, 4))
##' mcstate::array_flatten(x, 2:3)
##'
##' # array_flatten and array_reshape are each others' conceptual
##' # opposites:
##' y <- mcstate::array_flatten(x, 2:3)
##' identical(mcstate::array_reshape(y, 2, c(3, 4)), x)
array_flatten <- function(x, i) {
  dx <- dim(x)
  assert_integer(i)
  if (any(i < 1 | i > length(dx))) {
    stop(sprintf("Values of 'i' must be in [1, %d]", length(dx)))
  }
  if (length(i) < 2) {
    stop("i must be vector of at least length 2")
  }
  if (any(diff(i) != 1)) {
    stop("All values of 'i' must be consecutive integers")
  }
  dx[[i[[1L]]]] <- prod(dx[i])
  dx_new <- dx[seq_along(dx)[-i[-1]]]
  if (length(dx_new) == 1L) {
    dx_new <- NULL
  }
  dim(x) <- dx_new
  x
}


`array_last_dimension<-` <- function(x, i, value) { # nolint
  rank <- length(dim(x))
  if (rank == 2) {
    x[, i] <- value
  } else if (rank == 3) {
    x[, , i] <- value
  } else if (rank == 4) {
    x[, , , i] <- value
  } else {
    stop("Unexpected rank")
  }
  x
}


array_last_dimension <- function(x, i, drop = FALSE) {
  rank <- length(dim(x))
  if (rank == 2) {
    x[, i, drop = drop]
  } else if (rank == 3) {
    x[, , i, drop = drop]
  } else if (rank == 4) {
    x[, , , i, drop = drop]
  } else {
    stop("Unexpected rank")
  }
}


array_first_dimension <- function(x, i, drop = FALSE) {
  rank <- length(dim(x))
  if (rank == 2) {
    x[i, , drop = drop]
  } else if (rank == 3) {
    x[i, , , drop = drop]
  } else if (rank == 4) {
    x[i, , , , drop = drop]
  } else {
    stop("Unexpected rank")
  }
}


## This will be slower than above.
array_nth_dimension <- function(x, k, i, drop = FALSE) {
  rank <- length(dim(x))
  if (rank == 2) {
    expr <- quote(x[, , drop = FALSE])
  } else if (rank == 3) {
    expr <- quote(x[, , , drop = FALSE])
  } else if (rank == 4) {
    expr <- quote(x[, , , , drop = FALSE])
  } else {
    stop("Unexpected rank")
  }
  if (k < 1 || k > rank) {
    stop(sprintf("'k' must be in [1, %d]", rank))
  }

  expr[[k + 2]] <- quote(i)
  eval(expr)
}


array_from_list <- function(data, order = NULL) {
  len <- lengths(data)
  if (all(len == 0)) {
    return(NULL)
  }

  dimensions <- c(dim(data[[1]]) %||% length(data[[1]]), length(data))
  rank <- length(dimensions)

  stopifnot(
    length(unique(len)) == 1,
    is.null(order) || setequal(order, seq_len(rank)))

  vec <- unlist(data, FALSE, FALSE)
  if (rank == 2) {
    if (is.null(order) || order[[1]] == 1L) {
      return(matrix(vec, ncol = length(data)))
    } else {
      return(matrix(vec, nrow = length(data), byrow = TRUE))
    }
  }

  arr <- array(vec, dimensions)
  if (!is.null(order)) {
    arr <- aperm(arr, order)
  }

  arr
}
