mcstate_trajectories_discrete <- function(step, rate, state, predicted) {
  if (length(predicted) == 1L) {
    predicted <- rep(predicted, length(step))
  }
  ret <- list(step = step, rate = rate, state = state, predicted = predicted)
  class(ret) <- c("mcstate_trajectories_discrete", "mcstate_trajectories")
  ret
}


## There's some longer term tidying up this interface with the above,
## but the time issue is pretty fundamental unfortunately. Because the
## covid people depend on things as they are for now, I'm making a
## parallel class here.
mcstate_trajectories_continuous <- function(time, state, predicted) {
  if (length(predicted) == 1L) {
    predicted <- rep_len(predicted, length(step))
  }
  if (any(predicted)) {
    stop("predicted continuous trajectories not supported (mrc-3452, mrc-3453)")
  }
  ret <- list(time = time, state = state, predicted = predicted)
  class(ret) <- c("mcstate_trajectories_continuous", "mcstate_trajectories")
  ret
}


bind_mcstate_trajectories_discrete <- function(a, b) {
  stopifnot(inherits(a, "mcstate_trajectories_discrete"),
            inherits(b, "mcstate_trajectories_discrete"),
            last(a$step) == b$step[[1]],
            a$rate == b$rate,
            dim(a)[1:2] == dim(b)[1:2])

  step <- c(a$step, b$step[-1])
  if (length(dim(b$state)) == 3) {
    state <- array_bind(a$state, b$state[, , -1, drop = FALSE])
  } else {
    state <- array_bind(a$state, b$state[, , , -1, drop = FALSE])
  }
  rownames(state) <- rownames(b$state) %||% rownames(a$state)
  predicted <- c(a$predicted, b$predicted[-1])

  mcstate_trajectories_discrete(step, a$rate, state, predicted)
}


## Compatibility due to direct use in spimalot
mcstate_trajectories <- function(...) {
  .Deprecated("mcstate_trajectories_discrete")
  mcstate_trajectories_discrete(...)
}


bind_mcstate_trajectories <- function(...) {
  .Deprecated("bind_mcstate_trajectories")
  bind_mcstate_trajectories_discrete(...)
}
