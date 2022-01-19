`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}


vlapply <- function(x, fun, ...) {
  vapply(x, fun, logical(1), ...)
}


vnapply <- function(x, fun, ...) {
  vapply(x, fun, numeric(1), ...)
}


viapply <- function(x, fun, ...) {
  vapply(x, fun, integer(1), ...)
}


vcapply <- function(x, fun, ...) {
  vapply(x, fun, character(1), ...)
}


list_to_numeric <- function(x) {
  vnapply(x, identity)
}


data_frame <- function(...) {
  data.frame(..., stringsAsFactors = FALSE)
}


##' @importFrom stats rnorm
rmvnorm_generator <- function(vcv) {
  vcv <- unname(vcv)
  if (!isSymmetric(vcv, tol = sqrt(.Machine$double.eps))) {
    stop("vcv must be symmetric")
  }

  i_include <- !apply(vcv == 0, 1, all)

  if (!all(i_include)) {
    vcv <- vcv[i_include, i_include, drop = FALSE]
  }

  ev <- eigen(vcv, symmetric = TRUE)
  if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
    stop("vcv must be positive definite")
  }
  n <- nrow(vcv)
  res <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))

  if (all(i_include)) {
    function(mean, scale = 1.0) {
      mean + drop(rnorm(ncol(vcv)) %*% (res * sqrt(scale)))
    }
  } else {
    function(mean, scale = 1.0) {
      x <- numeric(length(i_include))
      x[i_include] <- drop(rnorm(ncol(vcv)) %*% (res * sqrt(scale)))
      mean + x
    }
  }
}


last <- function(x) {
  x[[length(x)]]
}


all_or_none <- function(x) {
  all(x) || !any(x)
}


squote <- function(x) {
  sprintf("'%s'", x)
}


r6_private <- function(x) {
  x[[".__enclos_env__"]]$private
}


set_into <- function(x, at, value) {
  x[at] <- value
  x
}


set_names <- function(x, nms) {
  names(x) <- nms
  x
}


str_collapse <- function(x) {
  paste0("{", paste0(x, collapse = ", "), "}")
}


recycle <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) == n) {
    x
  } else if (length(x) == 1L) {
    rep_len(x, n)
  } else {
    stop(sprintf("Invalid length for '%s', expected 1 or %d", name, n))
  }
}

is_3d_array <- function(x) {
  length(dim(x)) == 3
}


## Copied from ncol
nlayer <- function(x) {
  dim(x)[[3L]]
}


## Copied from NCOL
NLAYER <- function(x) { # nolint
  d <- dim(x)
  if (length(d) > 2L) {
    d[3L]
  } else {
    1L
  }
}


layernames <- function(x) {
  nms <- dimnames(x)
  if (!is.null(nms[[3L]])) {
    nms[[3L]]
  } else {
    NULL
  }
}


`layernames<-` <- function(x, value) { # nolint
  if (length(dim(x)) < 3) {
    stop("'x' has less than three dimensions")
  }

  nms <- dimnames(x)

  if (is.null(nms)) {
    if (is.null(value)) {
      stop("'value' cannot be NULL if 'dimnames(x)' is NULL")
    }
    nms <- vector("list", length(dim(x)))
  }

  if (is.null(value)) {
    nms[3L] <- list(NULL)
  } else {
    nms[[3L]] <- assert_scalar_character(value)
  }

  dimnames(x) <- nms
  x
}


set_layernames <- function(m, nms) {
  layernames(m) <- nms
  m
}


normalise <- function(x) {
  x / sum(x)
}


test_integer <- function(x, name = deparse(substitute(x)),
                         what = "integer") {
  if (!(is.integer(x))) {
    eps <- sqrt(.Machine$double.eps)
    usable_as_integer <- is.numeric(x) && (max(abs(round(x) - x)) < eps)
    if (!usable_as_integer) {
      return(FALSE)
    }
  }

  TRUE
}


## I've seen rare cases where the processx poll says that the process
## is ready, but the process is still alive (shutting down).  So we
## give a short grace window here.
callr_safe_result <- function(rs, grace = 2, dt = 0.1) {
  for (i in seq(0, grace, by = dt)) {
    if (!rs$is_alive()) {
      break
    }
    Sys.sleep(dt)
  }
  ## Failure here will throw and propagate nicely to the controlling
  ## process.
  rs$get_result()
}
