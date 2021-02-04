# mcstate 0.4.3

* Support for incrementally running a particle filter (up to some point in the time series) and forking these partial runs; see the `$begin_run` method on the particle filter (#78)

# mcstate 0.4.2

* Fix typo in `sir_models.Rmd`

# mcstate 0.4.1

* New `$fix()` method on `pmcmc_parameters` objects for fixing the value for a subset of parameters before running with `pmcmc` (#98)

# mcstate 0.4.0

* Compare functions no longer use (or accept) the `prev_state` argument and now use just the current model state. This requires that models compute things like "daily incidence" within model code but simplifies use with irregular time series (#94)

# mcstate 0.3.4

* Support for "compiled compare functions", introduced in `dust` 0.6.1 (#92)

# mcstate 0.3.1

* The particle filter can now return the entire model state at points during the run, with argument `save_restart` to `$run()` and method `$restart_state()` (#86)
* The `pmcmc` can returned sample restart state using the `save_restart` argument to `mcstate::pmcmc_control` which can be used to restart the pMCMC part way along the time series (see `vignette("restart")`)

# mcstate 0.3.0

* `pmcmc` is now controllable via a new `mcstate::pmcmc_control` object
* `pmcmc` can run chains in parallel using `callr`, by specifying `n_workers = n` for `n` greater than 1.

# mcstate 0.2.16

* `pmcmc` adds new `rerun_every` argument to rerun the particle filter unconditionally.
