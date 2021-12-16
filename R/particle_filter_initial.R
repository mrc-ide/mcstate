##' Create a suitable initial condition function from a set of restart
##' state.  This takes care of a few bookkeping and serialisation
##' details and returns a function appropriate to pass to
##' [mcstate::particle_filter] as `initial`.
##'
##' @title Create restart initial state
##'
##' @param state A matrix of state (rows are different states, columns
##'   are different realisations).  This is the form of a slice pulled
##'   from a restart.
##'
##' @return A function with arguments `info`, `n_particles` and `pars`
##'   that will sample, with replacement, a matrix of state suitable
##'   as a starting point for a particle filter.  The `info` and
##'   `pars` arguments are ignored.
##'
##' @export
particle_filter_initial <- function(state) {
  assert_is(state, "matrix")
  e <- new.env(parent = baseenv())
  e$state <- state
  with(e, function(info, n_particles, pars) {
    state[, sample.int(ncol(state), n_particles, replace = TRUE), drop = FALSE]
  })
}
