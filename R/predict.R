##' Run predictions from the results of [pmcmc()]. This
##' function can also be called by running [predict()] on
##' the object, using R's S3 dispatch.
##'
##' @title Run predictions from PMCMC
##'
##' @param object The results of running [pmcmc()] with
##'   `return_state = TRUE` (without this extra information,
##'   prediction is not possible)
##'
##' @param steps A vector of time steps to return predictions for. The
##'   first value must be the final value run in your simulation. An
##'   error will be thrown if you get this value wrong, look in
##'   `object$predict$step` (or the error message) for the
##'   correct value.
##'
##' @param n_threads The number of threads used in the simulation. If
##'   not given, we default to the value used in the particle filter
##'   that was used in the pmcmc.
##'
##' @param seed The random number seed (see [`particle_filter`]). The
##'   default value of `NULL` will seed the dust random number
##'   generator from R's random number generator. However, you can
##'   pick up from the same RNG stream used in the simulation if you
##'   pass in `seed = object$predict$seed`. However, do not do this if
##'   you are gong to run `pmcmc_predict()` multiple times the result
##'   will be identical. If you do want to call predict with this
##'   state multiple times you should call
##'   [dust::dust_rng_state_long_jump()] with a `times` argument of `i`
##'   for the `i`'th usage (i.e., once for the first usage,
##'   two times for the 2nd, etc).
##'
##' @param prepend_trajectories Prepend trajectories from the particle
##'   filter to the predictions created here.
##'
##' @export
pmcmc_predict <- function(object, steps, prepend_trajectories = FALSE,
                          n_threads = NULL, seed = NULL) {
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
  pars <- apply(object$pars, 1, object$predict$transform)
  index <- object$predict$index
  model <- object$predict$filter$model
  n_threads <- n_threads %||% object$predict$filter$n_threads

  ## NOTE: n_particles here is 1 because every particle gets a
  ## different state.
  mod <- model$new(pars, steps[[1]], 1L, n_threads = n_threads,
                   seed = seed, pars_multi = TRUE)
  mod$set_state(array_reshape(state, 2, c(1, length(pars))))
  mod$set_index(index)
  y <- mod$simulate(steps)
  y <- array_drop(y, 2L) # (state x particles(1) x pars x time)

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
