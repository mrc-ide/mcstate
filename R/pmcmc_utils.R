mcstate_pmcmc <- function(pars, probabilities, state, trajectories) {
  ret <- list(pars = pars,
              probabilities = probabilities,
              state = state,
              trajectories = trajectories)
  class(ret) <- "mcstate_pmcmc"
  ret
}


##' @export
format.mcstate_pmcmc <- function(x, ...) {
  if (is.null(x$state)) {
    str_state <- sprintf("  state: (not included)")
  } else {
    str_state <- sprintf("  state: %d x %d matrix of final states",
                         nrow(x$state), ncol(x$state))
  }
  if (is.null(x$trajectories)) {
    str_trajectories <- sprintf("  trajectories: (not included)")
  } else {
    str_trajectories <- sprintf(
      "  trajectories: %d x %d x %d array of particle trajectories",
      nrow(x$trajectories), ncol(x$trajectories), dim(x$trajectories)[[3]])
  }

  indent <- 4
  c(sprintf("<mcstate_pmcmc> (%d samples)", nrow(x$pars) - 1L),
    sprintf("  pars: %d x %d matrix of parameters", nrow(x$pars), ncol(x$pars)),
    strwrap(paste(colnames(x$pars), collapse = ", "),
            indent = indent, exdent = indent),
    sprintf("  probabilities: %d x %d matrix of log-probabilities",
            nrow(x$probabilities), ncol(x$probabilities)),
    strwrap(paste(colnames(x$probabilities), collapse = ", "),
            indent = indent, exdent = indent),
    str_state,
    str_trajectories)
}


##' @export
print.mcstate_pmcmc <- function(x, ...) {
  cat(paste0(format(x), "\n", collapse = ""))
  invisible(x)
}
