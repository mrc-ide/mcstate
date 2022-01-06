mcstate_pmcmc <- function(pars, probabilities, state, trajectories, restart,
                          predict, chain = NULL, iteration = NULL) {

  iteration <- iteration %||% seq.int(0, length.out = nrow(pars))

  nested <- length(dim(pars)) == 3

  ret <- list(nested = nested,
              chain = chain,
              iteration = iteration,
              pars = pars,
              probabilities = probabilities,
              state = state,
              trajectories = trajectories,
              restart = restart,
              predict = predict)
  class(ret) <- "mcstate_pmcmc"
  ret
}


##' @export
format.mcstate_pmcmc <- function(x, ...) {
  format_dims <- function(x) {
    paste(paste(dim(x), collapse = " x "),
          if (length(dim(x)) == 2) "matrix" else "array")
  }

  if (is.null(x$state)) {
    str_state <- sprintf("  state: (not included)")
  } else {
    str_state <- sprintf("  state: %s of final states", format_dims(x$state))
  }
  if (is.null(x$trajectories)) {
    str_trajectories <- sprintf("  trajectories: (not included)")
  } else {
    str_trajectories <- sprintf(
      "  trajectories: %s of particle trajectories",
      format_dims(x$trajectories$state))
  }

  if (is.null(x$restart)) {
    str_restart <- sprintf("  restart: (not included)")
  } else {
    str_restart <- sprintf(
      "  restart: %s of particle restart state",
      format_dims(x$restart$state))
  }

  if (is.null(x$chain)) {
    header <- sprintf("<mcstate_pmcmc> (%d samples)", nrow(x$pars))
  } else {
    header <- sprintf("<mcstate_pmcmc> (%d samples across %d chains)",
                      nrow(x$pars), length(unique(x$chain)))
  }

  indent <- 4
  if (isTRUE(x$nested)) { # isTRUE just for compatibility for now
    populations <- last(dimnames(x$pars))
    str_populations <- c(
      sprintf("  nested samples over %d populations:",
              length(populations)),
      strwrap(paste(populations, collapse = ", "),
              indent = indent, exdent = indent))
  } else {
    str_populations <- NULL
  }

  c(header,
    str_populations,
    sprintf("  pars: %s of parameters", format_dims(x$pars)),
    strwrap(paste(colnames(x$pars), collapse = ", "),
            indent = indent, exdent = indent),
    sprintf("  probabilities: %s of log-probabilities",
            format_dims(x$probabilities)),
    strwrap(paste(colnames(x$probabilities), collapse = ", "),
            indent = indent, exdent = indent),
    str_state,
    str_trajectories,
    str_restart)
}


##' @export
print.mcstate_pmcmc <- function(x, ...) {
  cat(paste0(format(x), "\n", collapse = ""))
  invisible(x)
}


## NOTE: we need to expose a 'force' argument here for testing, as
## otherwise under R CMD check the progress bar does not run.
pmcmc_progress <- function(n, show, force = FALSE) {
  if (show) {
    fmt <- "Step :current / :total [:bar] ETA :eta | :elapsedfull so far"
    t0 <- Sys.time()
    callback <- function(p) {
      message(sprintf("Finished %d steps in %s",
                      n, format(Sys.time() - t0, digits = 1)))
    }
    p <- progress::progress_bar$new(fmt, n, callback = callback, force = force)
    p$tick(0)
    p$tick
  } else {
    function() NULL
  }
}
