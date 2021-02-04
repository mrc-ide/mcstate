pmcmc_varied_parameter <- function(name, initial, min = -Inf, max = Inf,
                            discrete = FALSE, prior = NULL) {
  assert_scalar_character(name)
  assert_scalar_logical(discrete)
  assert_is(min, c("numeric", "integer"))
  assert_is(max, c("numeric", "integer"))

  lng <- length(initial)

  if (length(min) == 1) {
    min <- rep(min, lng)
  } else if (length(min) != lng) {
    stop(sprintf("Length of 'min' must be '1' or %s.", lng))
  }

  if (length(max) == 1) {
    max <- rep(max, lng)
  } else if (length(max) != lng) {
    stop(sprintf("Length of 'max' must be '1' or %s.", lng))
  }

  if (any(initial < min)) {
    stop(sprintf("'initial' must be >= 'min' (%s)", min))
  }
  if (any(initial > max)) {
    stop(sprintf("'initial' must be <= 'max' (%s)", max))
  }

  if (is.null(prior)) {
    prior <- function(p) { }
    body(prior) <- substitute(numeric(lng))
  } else {
    assert_is(prior, "function")
    value <- tryCatch(
      prior(initial),
      error = function(e)
        stop(sprintf(
          "Prior function for '%s' failed to evaluate initial value: %s",
          name, e$message)))
    if (!all(is.finite(value))) {
     stop(sprintf(
     "Prior function for '%s' returned a non-finite value on initial value(s)",
     name))
    }
    if (length(value) != lng) {
      stop(sprintf(
        "Prior function for '%s' should return %s values.",
        name, lng))
    }
  }

  ret <- list(name = name, initial = initial, min = min, max = max,
              discrete = discrete, prior = prior)
  class(ret) <- "pmcmc_varied_parameter"
  ret
}
