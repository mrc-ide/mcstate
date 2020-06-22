##' Run a grid search of the particle filter over two parameters
##'
##' @title Grid search over two parameters
##'
##' @param range \code{data.frame} describing the parameters to search.
##'  Must have columns name, min, max, n and target
##'
##' @param filter A \code{particle_filter} object to run
##'
##' @param n_particles Number of particles. Positive Integer. Default = 100
##'
##' @param tolerance Check around edges of the probability matrix
##'
##' @param run_params Passed to \code{particle_filter$run}
##'
##' @return List of beta and start date grid values, and
##'   normalised probabilities at each point
##'
##' @export
grid_search <- function(range, filter, n_particles, tolerance=1E-2,
                        run_params = NULL) {
  vars <- grid_search_validate_range(range)

  # Run the search, which returns a log-likelihood
  flat_log_ll <- vnapply(seq_len(nrow(vars$expanded)), function(i)
    filter$run2(n_particles, save_history = FALSE, index = vars$index,
                pars = vars$expanded[i, ], run_params = run_params))
  mat_log_ll <- matrix(
    flat_log_ll,
    nrow = length(vars$variables[[1]]),
    ncol = length(vars$variables[[2]]),
    byrow = FALSE
  )

  # Exponentiate elements and normalise to 1 to get probabilities
  prob_matrix <- exp(mat_log_ll)
  renorm_mat_ll <- prob_matrix / sum(prob_matrix)

  if (zero_boundary(renorm_mat_ll, tolerance = tolerance)) {
    message("Edges of the probability matrix not zero, check search range")
  }

  results <- list(
    vars = vars,
    x = vars$variables[1],
    y = vars$variables[2],
    mat_log_ll = mat_log_ll,
    renorm_mat_ll = renorm_mat_ll
  )

  class(results) <- "mcstate_scan"
  results
}


grid_search_validate_range <- function(range) {
  assert_is(range, "data.frame")
  msg <- setdiff(c("name", "min", "max", "n", "target"), names(range))
  if (length(msg) > 0L) {
    stop("Missing columns from 'range': ", paste(squote(msg), collapse = ", "))
  }

  if (nrow(range) != 2L) {
    stop("Expected exactly two rows in 'range'")
  }

  if (anyDuplicated(range$name)) {
    stop("Duplicate 'name' entries not allowed in 'range'")
  }

  targets <- c("step_start", "model_data", "pars_compare")
  err <- setdiff(range$target, targets)
  if (length(err) > 0L) {
    stop(sprintf("Invalid target %s: must be one of %s",
                 paste(squote(err), collapse = ", "),
                 paste(squote(targets), collapse = ", ")))
  }

  variables <- Map(seq, range$min, range$max, length.out = range$n)
  names(variables) <- range$name
  expanded <- do.call("expand.grid", variables)

  index <- lapply(targets, function(t) which(range$target == t))
  names(index) <- targets
  if (length(index$step_start) > 1L) {
    stop("At most one target may be 'step_start'")
  }

  list(range = range,
       variables = variables,
       expanded = expanded,
       index = index)
}

##' @export
plot.mcstate_scan <- function(x, ..., what = "likelihood", title = NULL) {
  if (what == "likelihood") {
    graphics::image(
      x = x$vars$variables[[1]], y = x$vars$variables[[2]], z = x$mat_log_ll,
      xlab = names(x$vars$variables)[1], ylab = names(x$vars$variables)[2],
      main = title)
  } else if (what == "probability") {
    graphics::image(
      x = x$vars$variables[[1]], y = x$vars$variables[[2]], z = x$renorm_mat_ll,
      xlab = names(x$vars$variables)[1], ylab = names(x$vars$variables)[2],
      main = title)
  }
}

# 2D only - checks edges of an array are less than a tolerance
zero_boundary <- function(array, tolerance) {
  zero_edges <-
    all(abs(array[1, ]) < tolerance) &&
    all(abs(array[, 1]) < tolerance) &&
    all(abs(array[nrow(array), ]) < tolerance) &&
    all(abs(array[, ncol(array)]) < tolerance)
  zero_edges
}
