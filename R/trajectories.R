mcstate_trajectories <- function(step, rate, state, predicted) {
  if (length(predicted) == 1L) {
    predicted <- rep(predicted, length(step))
  }
  ret <- list(step = step, rate = rate, state = state, predicted = predicted)
  class(ret) <- "mcstate_trajectories"
  ret
}


bind_mcstate_trajectories <- function(a, b) {
  stopifnot(inherits(a, "mcstate_trajectories"),
            inherits(b, "mcstate_trajectories"),
            last(a$step) == b$step[[1]],
            a$rate == b$rate,
            dim(a)[1:2] == dim(b)[1:2])

  step <- c(a$step, b$step[-1])
  state <- abind3(a$state, b$state[, , -1, drop = FALSE])
  rownames(state) <- rownames(b$state) %||% rownames(a$state)
  predicted <- c(a$predicted, b$predicted[-1])

  mcstate_trajectories(step, a$rate, state, predicted)
}
