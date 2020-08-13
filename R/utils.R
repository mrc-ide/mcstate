`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}


mcstate_file <- function(path) {
  system.file(path, package = "mcstate", mustWork = TRUE)
}


vlapply <- function(x, fun, ...) {
  vapply(x, fun, logical(1), ...)
}


vnapply <- function(x, fun, ...) {
  vapply(x, fun, numeric(1), ...)
}


list_to_numeric <- function(x) {
  vnapply(x, identity)
}


## Array-bind on 3rd dimension
abind3 <- function(a, b) {
  na <- dim(a)[3]
  nb <- dim(b)[3]
  nab <- dim(a)[1:2]
  ret <- array(NA_real_, c(nab, na + nb))
  ret[, , seq_len(na)] <- a
  ret[, , seq_len(nb) + na] <- b
  ret
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

  function(mean) {
    mean + drop(rnorm(ncol(vcv)) %*% res)
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
