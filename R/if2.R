##' @export
if2 <- function(pars, data, generator, compare, compare_pars, index,
                pars_sd, iterations, n_par_sets, cooling_target, # will become control
                n_threads = 1L, seed = NULL) {
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

  n_pars = length(pars)
  pars_vectors <- walk_initialise(pars, n_par_sets, n_steps, pars_sd)

  model <- generator$new(pars = pars_vectors$list_form,
                         step = steps[[1L]],
                         n_particles = 1L, n_threads = n_threads,
                         seed = seed, pars_multi = TRUE)
  if (!is.null(index)) {
    model$set_index(index()$state)
  }

  log_likelihood <- array(0, c(n_par_sets, iterations))
  if_pars <- array(NA_real_, c(c(n_par_sets, iterations, n_pars)))
  alpha_cool <- cooling_target ^ (1 / iterations)

  # TODO progress bar

  for (m in 1:iterations) {
    model$reset(pars = pars_vectors$list_form, steps[[1L]])
    for (t in seq(1L, n_steps)) {
      step_end <- steps[t, 2L]
      state <- model$run(step_end)[ , 1, , drop=TRUE]

      log_weights <- compare(state, data_split[[t]], compare_pars)

      if (!is.null(log_weights)) {
        weights <- scale_log_weights(log_weights)
        log_likelihood[, m] <- log_likelihood[, m] + weights$average
        if (log_likelihood == -Inf) {
          break
        }

        browser()
        kappa <- particle_resample(weights$weights)
        model$reorder(array(kappa, c(1, length(kappa))))
        pars_vectors <- pars_walk(pars_vectors$array_form[, kappa], pars_sd)
        model$set_pars(pars_vectors$list_form)
      }
    }
    sd <- sd * alpha_cool
    pars_final <- pars_vectors$array_form
    pars_vectors <- pars_walk(pars_final, pars_sd)
    if_pars[, m, ] <- pars_final
  }
  # outputs
  # pars: n_pars * iterations * n_par_sets
  # ll: iterations * n_par_sets
  list(log_likelihood = log_likelihood,
       if_pars = if_pars)
}

# TODO add plotting

# TODO add likelihood sample
# Run a particle filter at each point estimate at final state to get
# mean + standard error

pars_walk <- function(pars_vectors, pars_sd) {
  for (par_idx in 1:n_pars) {
    init <- pars_vectors[par_idx, ]
    walk <- rnorm(n_par_sets, init, pars_sd[[par_idx]])
    pars_vectors[par_idx, ] <- walk
  }
  list(array_form = pars_vectors, list_form = apply(pars_vectors, 2, list))
}

# Set up pars_series and set the first step
walk_initialise <- function(pars, n_par_sets, n_steps, pars_sd) {
  n_pars = length(pars)
  pars_vectors = array(NA_real_, c(n_pars, n_par_sets),
                       dimnames = list(names(pars), NULL))
  for (par_idx in 1:n_pars) {
    walk <- rep(pars[[par_idx]], n_pars)
    walk <- rnorm(n_par_sets, walk, pars_sd[[par_idx]])
    pars_vectors[par_idx, ] <- walk
  }
  list(array_form = pars_vectors, list_form = apply(pars_vectors, 2, list))
}