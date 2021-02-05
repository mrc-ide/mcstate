##' Control for [smc2]. This function constructs a list of options and
##' does some basic validation to ensure that the options will work
##' well together. Do not manually change the values in this
##' object. Do not refer to any argument except `n_parameter_sets` by
##' position as the order of the arguments may change in future.
##'
##' @title Control for SMC2
##'
##' @param n_parameter_sets The number of replicate parameter sets to
##'   simulate at once.
##'
##' @param degeneracy_threshold The degeneracy threshold. Once the
##'   effective sample size drops below `degeneracy_threshold *
##'   n_parameter_sets` the algorithm will rerun simulations from the
##'   beginning of the data and use these to replenish the particles.
##'
##' @param covariance_scaling A scaling factor to update variance
##'   covariance matrix of sampled parameters by
##'
##' @param progress Logical, indicating if a progress bar should be
##'   displayed, using [`progress::progress_bar`].
##'
##' @param save_trajectories Logical, indicating if particle
##'   trajectories should be saved during the simulation.
##'
##' @return A `smc2_control` object, which should not be modified once created.
##' @examples
##'
##' mcstate::smc2_control(100)
smc2_control <- function(n_parameter_sets, degeneracy_threshold = 0.5,
                         covariance_scaling = 0.5, progress = TRUE,
                         save_trajectories = FALSE) {
  ret <- list(n_parameter_sets = n_parameter_sets,
              degeneracy_threshold = degeneracy_threshold,
              covariance_scaling = covariance_scaling,
              save_trajectories = save_trajectories,
              progress = progress)
  class(ret) <- "smc2_control"
  ret
}
