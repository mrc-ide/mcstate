##' @export
if2 <- function(pars, data, generator, compare, compare_pars, index,
                pars_sd, iterations, n_par_sets, cooling_target, # will become control
                n_threads = 1L, seed = NULL, progress = TRUE) {
  # particle filter setup

  if (!is.null(index) && !is.function(index)) {
    stop("'index' must be function if not NULL")
  }

  assert_is(data, "particle_filter_data")
  assert_is(pars, "if2_parameters")

  # TODO check inputs

  data_split <- df_to_list_of_lists(data)

  steps <- unname(as.matrix(data[c("step_start", "step_end")]))
  n_steps <- nrow(steps)

  pars_matrix <- pars$walk_initialise(n_par_sets, pars_sd)
  n_pars <- nrow(pars_matrix)

  model <- generator$new(pars = pars$model(pars_matrix),
                         step = steps[[1L]],
                         n_particles = NULL, n_threads = n_threads,
                         seed = seed, pars_multi = TRUE)
  if (!is.null(index)) {
    model$set_index(index(mod$info())$run)
  }

  log_likelihood <- array(0, c(n_par_sets, iterations))
  if_pars <- array(NA_real_, c(n_pars, n_par_sets, iterations))
  alpha_cool <- cooling_target^(1 / iterations)

  p <- pmcmc_progress(iterations, progress)

  for (m in seq_len(iterations)) {
    p()
    model$reset(pars = pars$model(pars_matrix), steps[[1L]])
    ## TODO: refactor internal loop into function
    for (t in seq_len(n_steps)) {
      step_end <- steps[t, 2L]
      state <- model$run(step_end)

      log_weights <- compare(state, data_split[[t]], compare_pars)
      log_weights <- log_weights + pars$prior(pars_matrix)

      if (!is.null(log_weights)) {
        weights <- scale_log_weights(log_weights)
        log_likelihood[, m] <- log_likelihood[, m] + weights$average
        if (any(log_likelihood[, m] == -Inf)) {
          break
        }

        kappa <- particle_resample(weights$weights)
        model$reorder(kappa)
        pars_vectors <- pars$walk(pars_matrix[, kappa], pars_sd)
        model$set_pars(pars$model(pars_matrix))
      }
    }
    pars_sd <- pars_sd * alpha_cool
    pars_final <- pars_matrix
    pars_matrix <- pars$walk(pars_final, pars_sd)
    browser()
    if_pars[, , m] <- pars_final
  }
  # outputs
  # pars: n_pars * n_par_sets * iterations
  # ll: iterations * n_par_sets
  list(log_likelihood = log_likelihood, if_pars = if_pars)
}

# TODO add plotting

# TODO add likelihood sample
# Run a particle filter at each point estimate at final state to get
# mean + standard error

