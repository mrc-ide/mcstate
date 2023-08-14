# mcstate 0.9.13

* Particle filters now work with irregular and non-unit spaced time series data

# mcstate 0.9.11

* Continuous time (ODE) models can now use workers for running chains in parallel with `pmcmc`

# mcstate 0.9.9

* Basic inference now working with continuous time (ODE) models, via `particle_filter` and `pmcmc`

# mcstate 0.9.2

* Support for adaptive proposals for deterministic models

# mcstate 0.9.1

* Allow running a particle filter with multiple parameter sets and a single data set.
* The `nested` field on the particle filter class has been split into two logical fields: `has_multiple_parameters` and `has_multiple_data`

# mcstate 0.9.0

* Deprecated 'discrete' argument to parameters in favour of 'integer' - affects `if2_parameter`, `pmcmc_parameter`, `pmcmc_varied_parameter`, `smc2_parameter` 

# mcstate 0.8.4

* Compiled compare functions now supported in more places - `particle_deterministic` and multistage models (#177)

# mcstate 0.8.3

* Overhaul `mcstate::pmcmc_chains_*` to always use a file for communication, making them easier to understand and more robust (#179)
* New functions `mcstate::pmcmc_chains_cleanup` for removing files created by the above, and `mcstate::pmcmc_chains_collect` for automating collecting samples
* New, simpler, approach to pmcmc parallelisation which shares as much code with the above.

# mcstate 0.8.2

* Allow filtering of the pmcmc chains during running (dropping burnin and filtering) to reduce memory usage when collectin large trajectories
* pmcmc no longer retains the initial parameter values

# mcstate 0.8.1

* New argument to `mcstate::particle_filter` and `mcstate::particle_deterministic`, `constant_log_likelihood` which can be used to compute the probabilities of non-time series data (#185)

# mcstate 0.8.0

* Rework the "nested" support; this now returns output in a different dimension order.  Primarily this is an internal refactoring.
* Allow use of multistage parameters with deterministic models, and with nested parameters.
* Transform functions for multistage parameters now take `info` and not `model` as an argument, more in keeping with other functions.

# mcstate 0.7.3

* Allow multistage parameters to work with the "deterministic" particle
* New `mcstate::particle_deterministic_state` object for advanced use of the deterministic particle
* Deterministic particle loses the `run_many` method

# mcstate 0.7.2

* Multistage particle filters now cope with running data covering a subset of their stages
* Drop support for chnging initial step via particle filter initial function for deterministic and nested filters

# mcstate 0.7.1

* New helper function `mcstate::particle_filter_initial` for creating particle filter initial state functions from restart data.
* Drop support for changing initial step via the particle filter initial function

# mcstate 0.7.0

* Multi-stage particle filter implemented, allowing arbitrary changes to model structure during a particle filter run (#159)

# mcstate 0.6.13

* Allow saving restart from the deterministic filter (#153)

# mcstate 0.6.5

* Reduced overhead in parallel pmcmc with workers, and faster/less memory-hungry chain combination (#142)

# mcstate 0.6.4

* Allow the particle filter to terminate early if we would not be interested in the result. This is useful for `mcstate::pmcmc` which can use it to stop calculating a likelood that would be rejected. Primarily useful when running with relatively low numbers of particles and a high variance in the estimator (#138)

# mcstate 0.6.3

* Add support for running in "deterministic" mode with recent dust (#139)

# mcstate 0.6.0

* Add an iterated filtering method via `mcstate::if2` (#123)

# mcstate 0.5.13

* New functions `pmcmc_chains_prepare` and `pmcmc_chains_run` which can be used to manually schedule chains over different computing resourcess (#129)

# mcstate 0.5.12

* When `rerun_every` is specified, a new control parameter `rerun_control` can be used to make this stochastic rerun

# mcstate 0.5.11

* The particle filter can now run entirely in compiled code if supported by the model. This may give a small performance gain, particularly on very simple models, or of the model has an expensive compare function (#118)

# mcstate 0.5.9

* Add `nested_step_ratio` parameter to `pmcmc_control` for controlling the ratio of fixed:varied steps for nested pMCMC

# mcstate 0.5.5

* New array helper `mcstate::array_flatten` for unshaping an array

# mcstate 0.5.4

* Remove deprecated arguments to `pmcmc` (these were deprecated in 0.3.0) (#114)

# mcstate 0.5.3

* Bugfix in `index` for nested particle filters.

# mcstate 0.5.2

* Extend support of `pmcmc` to `pmcmc_parameters_nested` objects.

# mcstate 0.5.1

* Added `particle_filter_state_nested` and extended `particle_filter` to handle `pmcmc_parameters_nested` objects.

# mcstate 0.5.0

* Basic SMC^2 implementation (`smc2()`) as an alternative to pmcmc. This is very embryonic and the interface will change over future versions to support things like restarting and saving trajectories (#13)

# mcstate 0.4.8

* Bug fixes in `$proposal` method of `pmcmc_parameters_nested` for discrete and bounded parameters.

# mcstate 0.4.7

* Added helper methods `mcstate::array_bind`, `mcstate::array_reshape` and `mcstate::array_drop` to simplify some common array operations (#106)

# mcstate 0.4.6

* Added `pmcmc_varied_parameter` for parameters that can vary between different populations.
* Added `pmcmc_parameters_nested` to hold parameters that vary between populations (`pmcmc_varied_parameter`) and parameters that are the same (fixed) between populations (`pmcmc_parameter`).

# mcstate 0.4.4

* Fix performance regression added in 0.4.3

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
