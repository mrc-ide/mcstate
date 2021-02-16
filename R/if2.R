##' @export
if2 <- function(pars, data, generator, compare, compare_pars,
                pars_sd, iterations, n_par_sets, cooling_target, # will become control
                n_threads = 1L, seed = NULL) {
  # particle filter setup
  if (!is_dust_generator(model)) {
    stop("'model' must be a dust_generator")
  }
  if (!is.null(index) && !is.function(index)) {
    stop("'index' must be function if not NULL")
  }
  if (!is.null(initial) && !is.function(initial)) {
    stop("'initial' must be function if not NULL")
  }
  assert_is(data, "particle_filter_data")

  #TODO check inputs

  data_split <- df_to_list_of_lists(data)

  steps <- unname(as.matrix(data[c("step_start", "step_end")]))
  n_steps <- nrow(steps)
  pars_series <- pars_walk(pars, n_par_sets, steps, pars_sd)
  pars_list <- apply(pars_series[, , 1], 2, list)

  #TODO pars transforms
  model <- generator$new(pars = pars_list, step = steps[[1L]],
                         n_particles = 1L, n_threads = n_threads,
                         seed = seed, pars_multi = TRUE)
  if (!is.null(index)) {
    model$set_index(index)
  }

  alpha_cool <- cooling_target ^ (1 / iterations)

  for (m in 1:n_par_sets) {
    model$reset(pars = pars_list, steps[[1L]])
    log_likelihood <- 0
    for (t in seq(1L, n_steps)) {
      step_end <- steps[t, 2L]
      state <- model$run(step_end)

      log_weights <- compare(state, data_split[[t]], compare_pars)

      if (is.null(log_weights)) {
        log_likelihood_step <- NA_real_
      } else {
        weights <- scale_log_weights(log_weights)
        log_likelihood_step <- weights$average
        log_likelihood <- log_likelihood + log_likelihood_step
        if (log_likelihood == -Inf) {
          break
        }

        kappa <- particle_resample(weights$weights)
        model$reorder(kappa)
        # TODO this is inefficient as shuffles all series forwards
        # TODO better to generate steps here rather than at the start?
        for (t_forward in seq(t, n_steps)) {
          pars_series[, , t_forward] <- pars_series[, kappa, t_forward]
        }
        pars_list <- apply(pars_series[, , t + 1L], 2, list)
        model$set_pars(pars_list)
      }
    }
    sd <- sd * alpha_cool
    pars_series <- pars_walk(pars_list, n_par_sets, steps, pars_sd)
    pars_list <- apply(pars_series[, , 1], 2, list)
  }
  # TODO outputs
  # should be:
  # pars: n_pars * iterations * n_par_sets
  # ll: iterations * n_par_sets
}

pars_walk <- function(pars, n_par_sets, steps, sd) {
  n_pars <- length(pars)
  pars_series = array(NA_real_, c(n_pars, n_par_sets, steps),
                      dimnames = list(names(pars), NULL, NULL))
  for (par_idx in 1:n_pars) {
    walk <- rep(pars[[par_idx]], n_pars)
    for (step in 1:steps) {
      walk <- rnorm(n_par_sets, walk, sd[par_idx])
      pars_series[par_idx, , step] <- walk
    }
  }
  pars_series
}