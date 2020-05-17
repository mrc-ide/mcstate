`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}


mcstate_file <- function(path) {
  system.file(path, package = "mcstate", mustWork = TRUE)
}


vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}


drop_dim <- function(x, n) {
  dim(x) <- dim(x)[-n]
  x
}
