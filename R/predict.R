pmcmc_predict <- function(object, steps, n_threads = NULL, seed = 1L) {
  if (is.null(object$predict)) {
    stop("mcmc was run with return_state = FALSE, can't predict")
  }
  if (length(steps) < 2) {
    stop("At least two steps required for predict")
  }
  if (steps[[1]] != object$predict$step) {
    stop(sprintf("Expected steps[1] to be %d", object$predict$step))
  }

  state <- object$state
  data <- apply(object$pars, 1, object$predict$transform)
  index <- object$predict$index
  model <- object$predict$model
  n_threads <- n_threads %||% object$predict$n_threads

  dust::dust_simulate(model, steps, data, state, index, n_threads, seed)
}


##' @export
predict.mcstate_pmcmc <- function(object, steps, n_threads = NULL, seed = 1L,
                                  ...) {
  pmcmc_predict(object, steps, n_threads, seed)
}
