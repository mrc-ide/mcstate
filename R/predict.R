##' Run predictions from the results of [pmcmc()]. This function can
##' also be called by running [predict()] on the object, using R's S3
##' dispatch.
##'
##' @title Run predictions from PMCMC
##'
##' @param object The results of running [pmcmc()] with `return_state
##'   = TRUE` (without this extra information, prediction is not
##'   possible)
##'
##' @param steps A vector of time steps to return predictions for. The
##'   first value must be the final value run in your simulation. An
##'   error will be thrown if you get this value wrong, look in
##'   `object$predict$step` (or the error message) for the correct
##'   value.
##'
##' @param n_threads The number of threads used in the simulation. If
##'   not given, we default to the value used in the particle filter
##'   that was used in the pmcmc.
##'
##' @param seed The random number seed - a positive integer.
##'
##' @param prepend_trajectories Prepend trajectories from the particle
##'   filter to the predictions created here.
##'
##' @export
pmcmc_predict <- function(object, steps, prepend_trajectories = FALSE,
                          n_threads = NULL, seed = 1L) {
  if (is.null(object$predict)) {
    stop("mcmc was run with return_state = FALSE, can't predict")
  }
  if (length(steps) < 2) {
    stop("At least two steps required for predict")
  }
  if (steps[[1]] != object$predict$step) {
    stop(sprintf("Expected steps[1] to be %d", object$predict$step))
  }
  if (prepend_trajectories && is.null(object$trajectories)) {
    stop(paste("mcmc was run with return_trajectories = FALSE,",
               "can't prepend trajectories"))
  }

  state <- object$state
  data <- apply(object$pars, 1, object$predict$transform)
  index <- object$predict$index
  model <- object$predict$model
  n_threads <- n_threads %||% object$predict$n_threads

  y <- dust::dust_simulate(model, steps, data, state, index, n_threads, seed)

  res <- mcstate_trajectories(steps, object$predict$rate, y, TRUE)
  if (prepend_trajectories) {
    res <- bind_mcstate_trajectories(object$trajectories, res)
  }

  res
}


##' @export
predict.mcstate_pmcmc <- function(object, ...) {
  pmcmc_predict(object, ...)
}
