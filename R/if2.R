##' @export
if2 <- function(pars, data, generator, compare, compare_pars, index,
                pars_sd, iterations, n_par_sets, cooling_target, # will become control
                n_threads = 1L, seed = NULL, progress = TRUE) {
  # particle filter setup

  if (!is.null(index) && !is.function(index)) {
    stop("'index' must be function if not NULL")
  }

  assert_is(data, "particle_filter_data")

  # TODO check inputs

  # TODO model initial states

  # TODO parameter transforms

  data_split <- df_to_list_of_lists(data)

  steps <- unname(as.matrix(data[c("step_start", "step_end")]))
  n_steps <- nrow(steps)

  n_pars <- length(pars)
  pars_vectors <- walk_initialise(pars, n_par_sets, pars_sd)

  model <- generator$new(pars = pars_vectors$list,
                         step = steps[[1L]],
                         n_particles = NULL, n_threads = n_threads,
                         seed = seed, pars_multi = TRUE)
  if (!is.null(index)) {
    model$set_index(index(mod$info())$run)
  }

  log_likelihood <- array(0, c(n_par_sets, iterations))
  if_pars <- array(NA_real_, c(c(n_par_sets, iterations, n_pars)))
  alpha_cool <- cooling_target^(1 / iterations)

  p <- pmcmc_progress(iterations, progress)

  for (m in seq_len(iterations)) {
    p()
    ## TODO: Anything needed if initial time is set as a parameter?
    ## Feels unlikely to work.
    model$reset(pars = pars_vectors$list, steps[[1L]])
    ## TODO: refactor internal loop into function
    for (t in seq_len(n_steps)) {
      step_end <- steps[t, 2L]
      state <- model$run(step_end)

      log_weights <- compare(state, data_split[[t]], compare_pars)

      if (!is.null(log_weights)) {
        weights <- scale_log_weights(log_weights)
        log_likelihood[, m] <- log_likelihood[, m] + weights$average
        if (any(log_likelihood[, m] == -Inf)) {
          break
        }

        kappa <- particle_resample(weights$weights)
        model$reorder(kappa)
        pars_vectors <- pars_walk(pars_vectors$matrix[, kappa], pars_sd)
        model$set_pars(pars_vectors$list)
      }
    }
    pars_sd <- pars_sd * alpha_cool
    pars_final <- pars_vectors$matrix
    pars_vectors <- pars_walk(pars_final, pars_sd)
    if_pars[, m, ] <- pars_final
  }
  # outputs
  # pars: n_pars * iterations * n_par_sets
  # ll: iterations * n_par_sets
  list(log_likelihood = log_likelihood, if_pars = if_pars)
}

# TODO add plotting

# TODO add likelihood sample
# Run a particle filter at each point estimate at final state to get
# mean + standard error

pars_walk <- function(pars, pars_sd) {
  n_pars <- nrow(pars) # also length(pars_sd)
  n_par_sets <- ncol(pars)
  ## This could easily be done in one fell swoop by vectorising over
  ## pars/sd; we might want to explore this if only working on a
  ## subset here.
  for (par_idx in seq_len(n_pars)) {
    pars[par_idx, ] <- rnorm(n_par_sets, pars[par_idx, ], pars_sd[[par_idx]])
  }
  list(matrix = pars, list = apply(pars, 2, as.list))
}

# Set up pars_series and set the first step
walk_initialise <- function(pars_initial, n_par_sets, pars_sd) {
  n_pars <- length(pars_initial)
  pars <- matrix(NA_real_, n_pars, n_par_sets)
  rownames(pars) <- names(pars_initial)
  for (par_idx in seq_len(n_pars)) {
    pars[par_idx, ] <-
      rnorm(n_par_sets, pars_initial[[par_idx]], pars_sd[[par_idx]])
  }
  list(matrix = pars, list = apply(pars, 2, as.list))
}
