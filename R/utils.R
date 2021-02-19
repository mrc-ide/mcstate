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
  ev <- eigen(vcv, symmetric = TRUE)
  if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
    stop("vcv must be positive definite")
  }
  n <- nrow(vcv)
  res <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))

  function(mean, scale = 1.0) {
    mean + drop(rnorm(ncol(vcv)) %*% (res * sqrt(scale)))
  }
}


list_to_matrix <- function(data) {
  len <- lengths(data)
  stopifnot(all(len == len[[1]]))
  len <- len[[1L]]
  matrix(unlist(data, FALSE, FALSE), length(data), len, byrow = TRUE)
}


list_to_array <- function(data) {
  len <- lengths(data)
  stopifnot(all(len == len[[1L]]))
  array(unlist(data, FALSE, FALSE), c(dim(data[[1L]]), length(data)))
}


set_colnames <- function(m, nms) {
  colnames(m) <- nms
  m
}


last <- function(x) {
  x[[length(x)]]
}


df_to_list_of_lists <- function(x) {
  lapply(unname(split(x, seq_len(nrow(x)))), as.list)
}

groupeddf_to_list_of_lists <- function(x, group) {
  ## largely copied from dust::dust_data
  rows <- lapply(seq_len(nrow(x)), function(i) as.list(x[i, ]))
  group <- x[[group]]
  rows_grouped <- unname(split(rows, group))
  lapply(seq_len(nrow(x) / length(unique(group))),
         function(i) lapply(rows_grouped, "[[", i))
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

nlayer <- function(x) {
  dim(x)[[3L]]
}
