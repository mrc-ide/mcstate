mcstate_pmcmc <- function(iteration, pars, probabilities, state,
                          trajectories, restart, predict, chain = NULL) {
  nested <- length(dim(pars)) == 3

  ## So the option here would be to either store the full
  if (nrow(pars) == length(iteration)) {
    pars_index <- NULL
  } else if (is.null(chain)) {
    pars_index <- iteration
  } else {
    ## We make the simplifying assumption that we always include the
    ## last iteration, which is done for us.  That *won't* be true
    ## after filtering, but that drops the full parameters so that's
    ## ok.
    len <- unname(tapply(iteration, chain, max))
    stopifnot(nrow(pars) == sum(len))
    pars_index <- iteration + cumsum(c(0, len[-length(len)]))[chain]
  }

  ret <- list(nested = nested,
              chain = chain,
              iteration = iteration,
              pars = pars,
              probabilities = probabilities,
              state = state,
              trajectories = trajectories,
              restart = restart,
              predict = predict,
              pars_index = pars_index)
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


##' @export
`[[.mcstate_pmcmc` <- function(x, i, ...) { # nolint
  assert_scalar_character(i)
  if (i %in% c("pars", "probabilities")) {
    ret <- NextMethod("[[")
    index <- x$pars_index
    if (!is.null(index)) {
      ret <- array_first_dimension(ret, index)
    }
    ret
  } else if (i %in% c("pars_full", "probabilities_full")) {
    i <- sub("_full$", "", i)
    NextMethod("[[")
  } else {
    NextMethod("[[")
  }
}


##' @export
`$.mcstate_pmcmc` <- function(x, name) { # nolint
  x[[name]]
}


## NOTE: we need to expose a 'force' argument here for testing, as
## otherwise under R CMD check the progress bar does not run.
pmcmc_progress <- function(control, force = FALSE) {
  if (control$progress) {
    n_steps <- control$n_steps
    ## TODO: tidy up this check here!
    if (identical(control$progress_style, "noninteractive")) {
      p <- progress_percentage(control$n_steps)
      p(0)
      p
    } else {
      fmt <- "Step :current / :total [:bar] ETA :eta | :elapsedfull so far"
      t0 <- Sys.time()
      callback <- function(p) {
        message(sprintf("Finished %d steps in %s",
                        n_steps, format(Sys.time() - t0, digits = 1)))
      }
      p <- progress::progress_bar$new(fmt, n_steps, callback = callback,
                                      force = force)
      p$tick(0)
      p$tick
    }
  } else {
    function() NULL
  }
}


progress_percentage <- function(total) {
  force(total)
  i <- 0
  prev <- -Inf
  function(n = 1) {
    i <<- i + n
    p <- floor(i / total * 100) # avoid 0.5 issues, report on completed steps
    if (p != prev) {
      prev <<- p
      message(paste("progress:", i))
    }
  }
}
