##' Control for [mcstate::if2()]. This function constructs a list of
##' options and does some basic validation. Do not manually change the
##' values in this object. Do not refer to any argument by position as
##' the order of the arguments may change in future.
##'
##' @title Control for IF2
##'
##' @param pars_sd The initial standard deviation of parameter walks.
##'
##' @param iterations The number of IF2 iterations to run, across which
##'   the cooling is performed
##'
##' @param n_par_sets The number of parameter sets to walk (c.f. the
##'   population size)
##'
##' @param cooling_target A factor < 1 multiplying pars_sd, which will be
##'   reached by the end of the iterations, and approached geometrically
##'
##' @param progress Logical, indicating if a progress bar should be
##'   displayed, using [`progress::progress_bar`].
##'
##' @return An `if2_control` object, which should not be modified once
##'   created. Pass this into [mcstate::if2()]
##'
##' @export
##' @examples
##' mcstate::if2_control(list(beta = 0.2, gamma = 0.2), 100, 1000, 0.5)
if2_control <- function(pars_sd, iterations, n_par_sets, cooling_target,
                        progress = TRUE) {
  assert_list_of(pars_sd, "numeric")
  assert_named(pars_sd)
  iterations <- assert_scalar_positive_integer(iterations)
  n_par_sets <- assert_scalar_positive_integer(n_par_sets)
  assert_logical(progress)

  if (cooling_target >= 1 || cooling_target <= 0) {
    stop("'cooling_target' must be between 0 and 1 (non-inclusive)")
  }

  ret <- list(pars_sd = pars_sd,
              iterations = iterations,
              n_par_sets = n_par_sets,
              cooling_target = cooling_target,
              progress = progress)
  class(ret) <- "if2_control"
  ret
}
