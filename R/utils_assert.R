assert_is <- function(x, what, name = deparse(substitute(x))) {
  if (!inherits(x, what)) {
    stop(sprintf("'%s' must be a %s", name,
                 paste(what, collapse = " / ")), call. = FALSE)
  }
  invisible(x)
}


## Special version of above to cope with classed functions
assert_function <- function(x, name = deparse(substitute(x))) {
  if (!is.function(x)) {
    stop(sprintf("'%s' must be a function", name), call. = FALSE)
  }
  invisible(x)
}


assert_function_or_null <- function(x, name = deparse(substitute(x))) {
  if (!is.null(x) && !is.function(x)) {
    stop(sprintf("'%s' must be function if not NULL", name), call. = FALSE)
  }
  invisible(x)
}


assert_named <- function(x, unique = FALSE, name = deparse(substitute(x))) {
  if (is.null(names(x))) {
    stop(sprintf("'%s' must be named", name), call. = FALSE)
  }
  if (unique && any(duplicated(names(x)))) {
    stop(sprintf("'%s' must have unique names", name), call. = FALSE)
  }
}


assert_integer <- function(x, name = deparse(substitute(x)),
                           what = "integer") {
  if (!(is.integer(x))) {
    eps <- sqrt(.Machine$double.eps)
    usable_as_integer <- is.numeric(x) && (max(abs(round(x) - x)) < eps)
    if (!usable_as_integer) {
      stop(sprintf("'%s' must be an %s", name, what), call. = FALSE)
    }
    x <- as.integer(round(x))
  }
  invisible(x)
}


assert_logical <- function(x, name = deparse(substitute(x))) {
  if (!(is.logical(x))) {
    stop(sprintf("'%s' must be a logical", name), call. = FALSE)
  }
  invisible(x)
}


assert_character <- function(x, name = deparse(substitute(x))) {
  if (!(is.character(x))) {
    stop(sprintf("'%s' must be a character", name), call. = FALSE)
  }
  invisible(x)
}


assert_strictly_increasing <- function(x, name = deparse(substitute(x))) {
  if (any(diff(x) <= 0)) {
    stop(sprintf("'%s' must be strictly increasing", name), call. = FALSE)
  }
  invisible(x)
}


assert_scalar <- function(x, name = deparse(substitute(x))) {
  if (length(x) != 1L) {
    stop(sprintf("'%s' must be a scalar", name), call. = FALSE)
  }
  invisible(x)
}


assert_scalar_positive_integer <- function(x, allow_zero = FALSE,
                                           name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  x <- assert_integer(x, name)
  min <- if (allow_zero) 0 else 1
  if (x < min) {
    stop(sprintf("'%s' must be at least %d", name, min), call. = FALSE)
  }
  invisible(x)
}


assert_scalar_logical <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  assert_logical(x, name)
  invisible(x)
}


assert_scalar_character <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  assert_character(x, name)
  invisible(x)
}

assert_list_of <- function(x, class, name = deparse(substitute(x))) {
  if (!(is.list(x))) {
    stop(sprintf("'%s' must be a list", name), call. = FALSE)
  }
  if (!all(vlapply(x, inherits, what = class))) {
    stop(sprintf("Elements of '%s' must be in '%s'", name,
                 str_collapse(class)), call. = FALSE)
  }
  invisible(x)
}


assert_dimensions <- function(x, expected, name = deparse(substitute(x))) {
  dim_x <- dim(x) %||% length(x)
  if (length(dim_x) != length(expected) || !all(dim_x == expected)) {
    rank <- length(dim_x)
    type <- c("a vector", "a matrix", "an array")[rank]
    dim <- if (rank == 1) "length" else "dimensions"
    stop(sprintf(
      "Expected '%s' to be %s with %s %s",
      name, type, dim, paste(expected, collapse = " x ")))
  }
  invisible(x)
}


assert_dimnames <- function(x, expected, name = deparse(substitute(x))) {
  rank <- length(expected)
  dn_x <- if (rank == 1) list(names(x)) else dimnames(x)
  if (!is.null(dn_x) && !identical(dn_x, unname(expected))) {
    for (i in seq_along(expected)) {
      if (!is.null(dn_x[[i]]) && !identical(dn_x[[i]], expected[[i]])) {
        if (is.null(expected[[i]])) {
          if (rank == 1) {
            stop(sprintf("Expected '%s' to have no names", name))
          } else {
            stop(sprintf("Expected names of dimension %d of '%s' to be empty",
                         i, name))
          }
        } else {
          nms <- names(expected)
          values <- paste(squote(expected[[i]]), collapse = ", ")
          if (is.null(nms) || !nzchar(nms[[i]])) {
            target <- values
          } else {
            target <- sprintf("%s (%s)", nms[[i]], values)
          }
          if (rank == 1) {
            stop(sprintf("Expected names of '%s' to match %s", name, target))
          } else {
            stop(sprintf("Expected names of dimension %d of '%s' to match %s",
                         i, name, target))
          }
        }
      }
    }
  }
  if (rank == 1) {
    names(x) <- expected[[1]]
  } else {
    dimnames(x) <- unname(expected)
  }
  invisible(x)
}


match_value <- function(arg, choices, name = deparse(substitute(arg))) {
  assert_scalar_character(arg)
  if (!(arg %in% choices)) {
    stop(sprintf("'%s' must be one of %s",
                 name, paste(squote(choices), collapse = ", ")))
  }
  arg
}
