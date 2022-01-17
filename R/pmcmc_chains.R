pmcmc_chains_prepare <- function(path, pars, filter, control,
                                 initial = NULL) {
  path <- pmcmc_chains_path(path)
  assert_is(pars, c("pmcmc_parameters", "pmcmc_parameters_nested"))
  assert_is(filter, c("particle_filter", "particle_deterministic"))
  assert_is(control, "pmcmc_control")

  if (control$n_workers != 1) {
    stop("'n_workers' must be 1")
  }
  if (!control$use_parallel_seed) {
    stop("'use_parallel_seed' must be TRUE")
  }

  initial <- pmcmc_check_initial(initial, pars, control$n_chains)
  inputs <- filter$inputs()

  seed <- make_seeds(control$n_chains, inputs$seed, filter$model)
  dat <- list(pars = pars, initial = initial, filter = inputs,
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

  seed <- inputs$seed[[chain_id]]
  set.seed(seed$r) # likely problematic on CRAN...
  filter <- particle_filter_from_inputs(inputs$filter, seed$dust)
  filter$set_n_threads(control$n_threads)

  initial <- inputs$initial[[chain_id]]

  samples <- pmcmc_run_chain(inputs$pars, initial, filter, control)

  saveRDS(samples, path$results)

  path$results
}


pmcmc_chains_collect <- function(path) {
  control <- readRDS(pmcmc_chains_path(path)$control)
  path <- pmcmc_chains_path(path, seq_len(control$n_chains))

  ## TODO: better error message
  msg <- !file.exists(path$results)
  if (any(msg)) {
    stop(sprintf("Results missing for chains %s",
                 paste(which(msg), collapse = ", ")))
  }

  ## Simplest way first, highest memory use; later we will try and do
  ## something more clever here by building things up more nicely.
  samples <- lapply(path$results, readRDS)
  pmcmc_combine(samples = samples)
}


pmcmc_chains_cleanup <- function(path) {
  control <- readRDS(pmcmc_chains_path(path)$control)
  path <- pmcmc_chains_path(path, seq_len(control$n_chains))
  unlink(c(path$inputs, path$control, path$results))
  if (length(dir(path$root, all.files = TRUE, no.. = TRUE)) == 0) {
    unlink(path$root, recursive = TRUE)
  }
}


pmcmc_chains_path <- function(path, chain_id = NULL) {
  assert_scalar_character(path)
  list(root = path,
       inputs = file.path(path, "inputs.rds"),
       control = file.path(path, "control.rds"),
       results = file.path(path, sprintf("results_%d.rds", chain_id)))
}