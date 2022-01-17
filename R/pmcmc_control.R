##' Control for the pmcmc. This function constructs a list of options
##' and does some basic validation to ensure that the options will
##' work well together. Do not manually change the values in this
##' object. Do not refer to any argument except `n_steps` by position
##' as the order of the arguments may change in future.
##'
##' pMCMC is slow and you will want to parallelise it if you possibly
##' can. There are two ways of doing this which are discussed in some
##' detail in `vignette("parallelisation", package = "mcstate")`.
##'
##' @section Thinning the chain at generation:
##'
##' Generally it may be preferable to thin the chains after generation
##'   using [mcstate::pmcmc_thin] or [mcstate::pmcmc_sample].
##'   However, waiting that long can create memory consumption issues
##'   because the size of the trajectories can be very large.  To
##'   avoid this, you can thin the chains at generation - this will
##'   avoid creating large trajectory arrays, but will discard some
##'   information irretrivably.
##'
##' If either of the options `n_burnin` or `n_steps_retain` are provided,
##'   then we will subsample the chain at generation.
##'
##' * If `n_burnin` is provided, then the first `n_burnin` (of
##'   `n_steps`) samples is discarded.  This must be at most `n_steps`
##' * If `n_steps_retain` is provided, then we *evenly* sample out of
##'   the remaining samples.  The algorithm will try and generate a
##'   sensible set here, and will always include the last sample of
##'   `n_steps` but may not always include the first post-burnin
##'   sample.  An error will be thrown if a suitable sampling is not
##'   possible (e.g., if `n_steps_retain` is larger than `n_steps -
##'   n_burnin`
##'
##' If either of `n_burnin` or `n_steps_retain` is provided, the
##'   resulting samples object will include the full set of parameters
##'   and probabilities sampled, along with an index showing how they
##'   relate to the filtered samples.
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
##' @param n_steps_each If using workers (i.e., `n_workers > 1`), the
##'   number of steps to run in each "chunk" on each worker before
##'   reporting back to the main process. Increasing this will make
##'   progress reporting less frequent and reduce some communication
##'   overhead (though the overhead is likely to be trivial in any
##'   real application). Decreasing this will give more frequent
##'   process reporting and if `n_threads_total` is given will allow
##'   for more rapid re-allocation of unused cores once chains start
##'   finishing. The default, if not given and if `n_workers > 1` is
##'   to use 10% of `n_steps`.
##'
##' @param n_threads_total The total number of threads (i.e., cores)
##'   the total number of threads/cores to use. If `n_workers` is
##'   greater than 1 then these threads will be divided evenly across
##'   your workers at first and so `n_threads_total` must be an even
##'   multiple of `n_workers`. If chains finish at different times
##'   (including if `n_chains` is not a multiple of `n_workers`) then
##'   these threads/cores will be reallocated across workers that are
##'   still going. If `n_workers` is 1 (i.e., running in parallel) and
##'   `n_threads_total` is not given (i.e., `NULL`) we will use the
##'   number of threads specified in the particle filter
##'   creation. Otherwise this value overrides the value in the
##'   particle filter.
##'
##' @param n_workers Number of "worker" processes to use to run chains
##'   in parallel. This must be at most `n_chains` and is recommended
##'   to be a divisor of `n_chains`. If `n_workers` is 1, then chains
##'   are run in series (i.e., one chain after the other). See the
##'   parallel vignette (`vignette("parallelisation", package =
##'   "mcstate")`) for more details about this approach.
##'
##' @param rerun_every Optional integer giving the frequency at which
##'   we should rerun the particle filter on the current "accepted"
##'   state.  The default for this (`Inf`) will never rerun this
##'   point, but if you set to 100, then every 100 steps we run the
##'   particle filter on both the proposed *and* previously accepted
##'   point before doing the comparison.  This may help "unstick"
##'   chains, at the cost of some bias in the results.
##'
##' @param rerun_random Logical, controlling the behaviour of
##'   rerunning (when `rerun_every` is finite). The default value of
##'   `FALSE` will rerun the filter deterministically at a fixed
##'   number of iterations (given by `rerun_every`). If `TRUE`, then
##'   we stochastically rerun each step with probability of `1 /
##'   rerun_every`. This gives the same expected number of MCMC steps
##'   between reruns but a different pattern.
##'
##' @param use_parallel_seed Logical, indicating if seeds should be
##'   configured in the same way as when running workers in parallel
##'   (with `n_workers > 1`).  Set this to `TRUE` to ensure
##'   reproducibility if you use this option sometimes (but not
##'   always). This option only has an effect if `n_workers` is 1.
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
##' @param save_restart An integer vector of time points to save
##'   restart information for; this is in addition to `save_state`
##'   (which saves the final model state) and saves the full model
##'   state.  It will use the same trajectory as `save_state` and
##'   `save_trajectories`. Note that if you use this option you will
##'   end up with lots of model states and will need to process them
##'   in order to actually restart the pmcmc or the particle filter
##'   from this state. The integers correspond to the *time* variable
##'   in your filter (see [mcstate::particle_filter] for more
##'   information).
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
##' @param nested_step_ratio Either integer or 1/integer, which specifies the
##'   ratio of fixed:varied steps in a nested pMCMC. For example `3` would run
##'   3 steps proposing fixed parameters only and then 1 step proposing varied
##'   parameters only; whereas `1/3` would run 3 varied steps
##'   for every 1 fixed step. The default value of `1` runs an equal number of
##'   iterations updating the fixed and varied parameters. Sensible choices
##'   of this parameter may depend on the true ratio of fixed:varied parameters
##'   or on desired run-time, for example updating fixed parameters is
##'   quicker so more varied steps could be more efficient.
##'
##' @param nested_update_both If `FALSE` (default) then alternates
##'   between proposing fixed and varied parameter updates according
##'   to the ratio in `nested_step_ratio`. If `TRUE` then proposes
##'   fixed and varied parameters simultaneously and collectively
##'   accepts/rejects them, `nested_step_ratio` is ignored.
##'
##' @param filter_early_exit Logical, indicating if we should allow
##'   the particle filter to exit early for points that will not be
##'   accepted. Only use this if your log-likelihood never increases
##'   between steps. This will the the case where your likelihood
##'   calculation is a sum of discrete normalised probability
##'   distributions, but may not be for continuous distributions!
##'
##' @param n_burnin Optionally, the number of points to discard as
##'   burnin.  This happens separately to the burnin in
##'   [mcstate::pmcmc_thin] or [mcstate::pmcmc_sample].  See Details.
##'
##' @param n_steps_retain Optionally, the number of samples to retain from
##'   the `n_steps - n_burnin` steps.  See Details.
##'
##' @return A `pmcmc_control` object, which should not be modified
##'   once created.
##'
##' @export
##'
##' @examples
##' mcstate::pmcmc_control(1000)
##'
##' # Suppose we have a fairly large node with 16 cores and we want to
##' # run 8 chains. We can use all cores for a single chain and run
##' # the chains sequentially like this:
##' mcstate::pmcmc_control(1000, n_chains = 8, n_threads_total = 16)
##'
##' # However, on some platforms (e.g., Windows) this may only realise
##' # a 50% total CPU use, in which case you might benefit from
##' # splitting these chains over different worker processes (2-4
##' # workers is likely the largest useful number).
##' mcstate::pmcmc_control(1000, n_chains = 8, n_threads_total = 16,
##'                        n_workers = 4)
pmcmc_control <- function(n_steps, n_chains = 1L, n_threads_total = NULL,
                          n_workers = 1L, n_steps_each = NULL,
                          rerun_every = Inf, rerun_random = FALSE,
                          use_parallel_seed = FALSE,
                          save_state = TRUE, save_restart = NULL,
                          save_trajectories = FALSE, progress = FALSE,
                          nested_step_ratio = 1, nested_update_both = FALSE,
                          filter_early_exit = FALSE,
                          n_burnin = NULL, n_steps_retain = NULL) {
  assert_scalar_positive_integer(n_steps)
  assert_scalar_positive_integer(n_chains)
  assert_scalar_positive_integer(n_workers)

  ## Leave this be for now
  n_steps_each <- n_steps
  assert_scalar_positive_integer(n_steps_each)

  if (is.null(n_threads_total)) {
    n_threads <- 1L
    n_threads_total <- n_workers
  } else {
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
    n_threads <- n_threads_total / n_workers
  }

  if (!identical(unname(rerun_every), Inf)) {
    assert_scalar_positive_integer(rerun_every)
  }
  assert_scalar_logical(rerun_random)

  assert_scalar_logical(use_parallel_seed)
  assert_scalar_logical(save_state)
  assert_scalar_logical(save_trajectories)
  assert_scalar_logical(progress)
  assert_scalar_logical(filter_early_exit)

  if (n_chains < n_workers) {
    stop(sprintf("'n_chains' (%d) is less than 'n_workers' (%d)",
                 n_chains, n_workers))
  }

  if (n_chains %% n_workers != 0) {
    ## TODO: we should improve this message
    message("Your last workers will be underpowered")
  }

  if (!is.null(save_restart)) {
    ## possibly assert_integer(save_restart)?
    assert_strictly_increasing(save_restart)
  }

  assert_scalar_logical(nested_update_both)
  ok <- test_integer(nested_step_ratio) || test_integer(1 / nested_step_ratio)
  if (!ok) {
    stop(sprintf("Either 'nested_step_ratio' (%g) or 1/'nested_step_ratio'
                          must be an integer", nested_step_ratio))
  }

  filter <- pmcmc_filter_on_generation(n_steps, n_burnin, n_steps_retain)

  ret <- list(n_steps = n_steps,
              n_chains = n_chains,
              n_workers = n_workers,
              n_steps_each = n_steps_each,
              n_threads = n_threads,
              n_threads_total = n_threads_total,
              rerun_every = rerun_every,
              rerun_random = rerun_random,
              use_parallel_seed = use_parallel_seed,
              save_state = save_state,
              save_restart = save_restart,
              save_trajectories = save_trajectories,
              progress = progress,
              filter_early_exit = filter_early_exit,
              nested_update_both = nested_update_both,
              nested_step_ratio = nested_step_ratio)
  ret[names(filter)] <- filter

  class(ret) <- "pmcmc_control"
  ret
}


## What do we do here about our starting point?  Probably best to just
## forget that for now really, as it's not really a sample...
pmcmc_filter_on_generation <- function(n_steps, n_burnin, n_steps_retain) {
  n_burnin <- assert_scalar_positive_integer(n_burnin %||% 0, TRUE)
  if (n_burnin >= n_steps) {
    stop("'n_burnin' cannot be greater than or equal to 'n_steps'")
  }
  n_steps_possible <- n_steps - n_burnin
  n_steps_retain <- assert_scalar_positive_integer(
    n_steps_retain %||% n_steps_possible)
  if (n_steps_retain > n_steps_possible) {
    stop(sprintf(
      "'n_steps_retain' is too large, max possible is %d but given %d",
      n_steps_possible, n_steps_retain))
  }

  ## Now, compute the step ratio:
  n_steps_every <- floor(n_steps_possible / n_steps_retain)
  seq(to = n_steps, length.out = n_steps_retain, by = n_steps_every)

  if (n_steps_every == 1) {
    ## If we've dropped more than 5% of the chain this probably means
    ## that they're out of whack.  Need a good explanation here
    ## though.
    if ((n_steps_possible - n_steps_retain) / n_steps_possible > 0.05) {
      stop(paste("'n_steps_retain' is too large to skip any samples, and",
                 "would result in just increasing 'n_burnin' by more than",
                 "5% of your post-burnin samples. Please adjust 'n_steps'",
                 "'n_burnin' or 'n_steps_retain' to make your intentions",
                 "clearer"))
    }
  }

  ## Back calculate the actual number of burnin steps to take:
  n_burnin <- n_steps - n_steps_every * (n_steps_retain - 1) - 1

  list(n_burnin = n_burnin,
       n_steps_retain = n_steps_retain,
       n_steps_every = n_steps_every)
}


pmcmc_check_control <- function(control) {
  ## An error here would mean history saving would fail in peculiar ways
  ok <- control$n_steps ==
    control$n_burnin + (control$n_steps_retain - 1) * control$n_steps_every + 1
  if (!ok) {
    stop(paste("Corrupt pmcmc_control (n_steps/n_steps_retain/n_burnin),",
               "perhaps you modified it after creation?"))
  }
}
