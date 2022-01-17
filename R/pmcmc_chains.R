pmcmc_chains_prepare <- function(path, pars, filter, control,
                                 initial = NULL) {
  path <- pmcmc_chains_path(path)
  assert_is(pars, c("pmcmc_parameters", "pmcmc_parameters_nested"))
  assert_is(filter, "particle_filter")
  assert_is(control, "pmcmc_control")

  ## These feel like we might modify them?
  ## if (control$n_workers != 1) {
  ##   stop("'n_workers' must be 1")
  ## }
  ## if (!control$use_parallel_seed) {
  ##   stop("'use_parallel_seed' must be TRUE")
  ## }

  initial <- pmcmc_check_initial(initial, pars, control$n_chains)

  seed <- make_seeds(control$n_chains, filter$inputs()$seed, filter$model)

  dat <- list(pars = pars, initial = initial, filter = filter,
              control = control, seed = seed)
  class(dat) <- "pmcmc_inputs"

  dir.create(path$root, FALSE, TRUE)
  ## TODO: avoids warnings about packages not being available on load,
  ## but we might erase those differently....
  suppressWarnings(saveRDS(dat, path$inputs))
  saveRDS(control, path$control)

  path$root
}


pmcmc_chains_run <- function(chain_id, path) {
  assert_scalar_positive_integer(chain_id)
  path <- pmcmc_chains_path(path, chain_id)

  inputs <- readRDS(path$inputs)
  assert_is(inputs, "pmcmc_inputs")

  control <- inputs$control

  if (chain_id < 1 || chain_id > control$n_chains) {
    stop(sprintf("'chain_id' must be an integer in 1..%d",
                 control$n_chains))
  }

  samples <- pmcmc_run_chain(chain_id, inputs$pars, inputs$initial,
                             inputs$filter, control, inputs$seed)

  saveRDS(samples, path$results)

  path$results
}


pmcmc_chains_collect <- function(path) {
  control <- readRDS(pmcmc_chains_path(path)$control)
  path <- pmcmc_chains_path(path, seq_len(control$n_chains))

  ## TODO: better error message
  if (!all(file.exists(path$results))) {
    stop("Some results missing")
  }

  ## Simplest way first, highest memory use; later we will try and do
  ## something more clever here by building things up more nicely.
  samples <- lapply(path$results, readRDS)
  pmcmc_combine(samples = samples)
}


pmcmc_chains_path <- function(path, chain_id = NULL) {
  assert_scalar_character(path)
  list(root = path,
       inputs = file.path(path, "inputs.rds"),
       control = file.path(path, "control.rds"),
       results = file.path(path, sprintf("results_%d.rds", chain_id)))
}
