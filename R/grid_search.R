##' Run a grid search of the particle filter over two parameters
##'
##' @title Grid search over two parameters
##'
##' @param range \code{data.frame} describing the parameters to search.
##'  Must have columns name, min, max, n and target
##'
##' @param filter A \code{particle_filter} object to run
##'
##' @param tolerance Check around edges of the probability matrix
##'
##' @return List of beta and start date grid values, and
##'   normalised probabilities at each point
##'
##' @export
grid_search <- function(range, filter, tolerance = 0.01) {
  vars <- grid_search_validate_range(range)

  # Run the search, which returns a log-likelihood
  flat_log_ll <- vnapply(seq_len(nrow(vars$expanded)), function(i) {
    pars <- vars$expanded[i, ]
    names(pars) <- colnames(vars$expanded)
    filter$run2(pars, vars$index)
  })

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
  # Defined in pmcmc.R
  index <- validate_range(range,
    c("name", "min", "max", "n", "target"))

  if (nrow(range) != 2L) {
    stop("Expected exactly two rows in 'range'")
  }

  variables <- Map(seq, range$min, range$max, length.out = range$n)
  names(variables) <- range$name
  expanded <- do.call("expand.grid", variables)

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
