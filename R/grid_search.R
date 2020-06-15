grid_search <- function(state, range, filter, n_particles) {
  vars <- grid_search_validate_range(range)

  res <- vnapply(seq_len(nrow(vars$expanded)), function(i)
    filter$run2(state, n_particles, vars$expanded[i, ], vars$index))
}


grid_search_validate_range <- function(range) {
  assert_is(range, "data.frame")
  msg <- setdiff(c("name", "min", "max", "n", "target"), names(range))
  if (length(msg) > 0L) {
    stop("Missing columns from 'range': ", paste(squote(msg), collapse = ", "))
  }
  targets <- c("step_start", "pars_model", "pars_compare")
  err <- setdiff(range$target, targets)
  if (length(err) > 0L) {
    stop(sprintf("Invalid target %s: must be one of %s",
                 paste(squote(err), collapse = ", "),
                 paste(squote(target), collapse = ", ")))
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
