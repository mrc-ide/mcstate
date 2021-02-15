array_bind <- function(..., arrays = list(...)) {
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
  d <- matrix(unlist(d), r)

  d_shared <- d[-r , , drop = FALSE]
  d_diff <- d_shared - d_shared[, 1] != 0
  if (any(d_diff)) {
    err <- which(d_diff, TRUE)
    err <- vcapply(split(unname(err[, 1]), err[, 2]), paste, collapse = ", ")
    detail <- paste(sprintf("array %s (dimension %s)", names(err), unname(err)),
                 collapse = ", ")
    stop(paste("Incompatible dimension arrays:", detail))
  }

  dn <- dimnames(arrays[[1L]])
  if (!is.null(dn)) {
    dn_last <- lapply(arrays, function(x) dimnames(x)[[r]])
    if (all(lengths(dn_last) > 0)) {
      dn[[r]] <- unlist(dn_last)
    } else {
      dn[r] <- list(NULL)
    }
  }

  array(unlist(arrays), c(d[-r, 1], sum(d[r, ])), dimnames = dn)
}


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

  ## The actual reshape is easy:
  dim(x) <- append(dx[-i], d, i - 1L)
  x
}
