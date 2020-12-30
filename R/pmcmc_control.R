##' Control for the pmcmc. This function constructs a list of options
##' and does some basic validation to ensure that the options will
##' work well together. Do not manually change the values in this
##' object. Do not refer to any argument except `n_steps` by position
##' as the order of the arguments may change in future.
##'
##' @title Control for the pmcmc
##'
##' @param n_steps Number of MCMC steps to run. This is the only
##'   required argument.
##'
##' @param n_chains Optional integer, indicating the number of chains
##'   to run. If more than one then we run a series of chains and
##'   merge them with [pmcmc_combine()]. Chains are run in series,
##'   with the same filter if `n_workers` is 1, or run in parallel
##'   otherwise.
##'
##' @param n_workers Number of "worker" processes to use to run chains
##'   in parallel. This must be greater than `n_chains` and is
##'   recommended to be a divisor of `n_chains`. If `n_workers` is 1,
##'   then chains are run in series (i.e., one chain after the other).
##'
##' @param n_steps_each If using workers (i.e., `n_workers > 1`), the
##'   number of steps to run in each "chunk" on each worker before
##'   reporting back to the main process. Increasing this will make
##'   progress reporting less frequent and reduce some communication
##'   overhead (though the overhead is likely to be trivial in any
##'   real application). Decreasing this will give more frequent
##'   process reporting and if `n_threads_total` is given will allow
##'   for more rapid re-allocation of unused cores once chains start
##'   finishing. The default, if not given and if `n_workers > 1` is
##'   to use `n_steps` which gives *no* within-chain feedback and no
##'   option to reallocate cores.
##'
##' @param n_threads_total If using workers (i.e., `n_workers > 1`),
##'   the total number of threads/cores to use. These threads will be
##'   divided evenly across your workers at first and so
##'   `n_threads_total` must be an even multiple of `n_workers`. If
##'   chains finish at different times (including if `n_chains` is not
##'   a multiple of `n_workers`) then these threads/cores will be
##'   reallocated across workers, provided that `n_steps_each` is
##'   given.
##'
##' @param rerun_every Optional integer giving the frequency at which
##'   we should rerun the particle filter on the current "accepted"
##'   state.  The default for this (`Inf`) will never rerun this
##'   point, but if you set to 100, then every 100 steps we run the
##'   particle filter on both the proposed *and* previously accepted
##'   point before doing the comparison.  This may help "unstick"
##'   chains, at the cost of some bias in the results.
##'
##' @param save_state Logical, indicating if the state should be saved
##'   at the end of the simulation. If `TRUE`, then a single
##'   randomly selected particle's state will be collected at the end
##'   of each MCMC step. This is the full state (i.e., unaffected by
##'   and `index` used in the particle filter) so that the
##'   process may be restarted from this point for projections.  If
##'   `save_trajectories` is `TRUE` the same particle will
##'   be selected for each. The default is `TRUE`, which will
##'   cause `n_state` * `n_steps` of data to be output
##'   alongside your results. Set this argument to `FALSE` to
##'   save space, or use [pmcmc_thin()] after running the
##'   MCMC.
##'
##' @param save_trajectories Logical, indicating if the particle
##'   trajectories should be saved during the simulation. If `TRUE`,
##'   then a single randomly selected particle's trajectory will be
##'   collected at the end of each MCMC step.  This is the filtered
##'   state (i.e., using the `state` component of `index` provided to
##'   the particle filter).  If `save_state` is `TRUE` the same
##'   particle will be selected for each.
##'
##' @param progress Logical, indicating if a progress bar should be
##'   displayed, using [`progress::progress_bar`].
##'
##' @return A `pmcmc_control` object, which should not be modified
##'   once created.
##'
##' @export
##'
##' @examples
##' mcstate::pmcmc_control(10)
pmcmc_control <- function(n_steps, n_chains = 1L, n_workers = 1L,
                          n_steps_each = NULL, n_threads_total = NULL,
                          rerun_every = Inf, save_state = TRUE,
                          save_trajectories = FALSE, progress = FALSE) {
  assert_scalar_positive_integer(n_steps)
  assert_scalar_positive_integer(n_chains)
  assert_scalar_positive_integer(n_workers)
  if (is.null(n_steps_each) || n_workers == 1L) {
    n_steps_each <- n_steps
  } else {
    assert_scalar_positive_integer(n_steps_each)
  }
  if (!is.null(n_threads_total)) {
    assert_scalar_positive_integer(n_threads_total)
    if (n_threads_total < n_workers) {
      stop(sprintf("'n_threads_total' (%d) is less than 'n_workers' (%d)",
                   n_threads_total, n_workers))
    }
    if (n_threads_total %% n_workers != 0) {
      stop(sprintf(
        "'n_threads_total' (%d) is not a multiple of 'n_workers' (%d)",
        n_threads_total, n_workers))
    }
  }

  if (!identical(unname(rerun_every), Inf)) {
    assert_scalar_positive_integer(rerun_every)

  }

  assert_scalar_logical(save_state)
  assert_scalar_logical(save_trajectories)
  assert_scalar_logical(progress)

  if (n_chains < n_workers) {
    stop(sprintf("'n_chains' (%d) is less than 'n_workers' (%d)",
                 n_chains, n_workers))
  }

  ret <- list(n_steps = n_steps,
              n_chains = n_chains,
              n_workers = n_workers,
              n_steps_each = n_steps_each,
              n_threads_total = n_threads_total,
              rerun_every = rerun_every,
              save_state = save_state,
              save_trajectories = save_trajectories,
              progress = progress)
  class(ret) <- "pmcmc_control"
  ret
}
