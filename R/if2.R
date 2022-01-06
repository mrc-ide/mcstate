##' Create an IF2 object for running and interacting with an IF2
##'   inference.
##'
##' See:
##' Ionides EL, Nguyen D, Atchadé Y, Stoev S, King AA (2015).
##' "Inference for Dynamic and Latent Variable Models via Iterated,
##' Perturbed Bayes Maps."
##' PNAS, 112(3), 719–724. https://doi.org/10.1073/pnas.1410597112.
##'
##' @title Run iterated filtering (IF2 algorithm)
##'
##' @param pars An [mcstate::if2_parameters] object, describing the
##'   parameters that will be varied in the simulation, and the method
##'   of transformation into model parameters.
##'
##' @param filter A [mcstate::particle_filter] object. We don't use
##'   the particle filter directly (except for sampling with
##'   `mcstate::if2_sample`) but this shares so much validation that
##'   it's convenient.  Be sure to set things like the seed and number
##'   of threads here if you want to use anything other than the
##'   default.
##'
##' @param control An [mcstate::if2_control()] object
##'
##' @return An object of class `if2_fit`, which contains the sampled
##'   parameters (over time) and their log-likelihoods
##'
##' @export
##' @importFrom R6 R6Class
##' @examples
##' # A basic SIR model used in the particle filter example
##' gen <- dust::dust_example("sir")
##'
##' # Some data that we will fit to, using 1 particle:
##' sir <- gen$new(pars = list(), step = 0, n_particles = 1)
##' dt <- 1 / 4
##' day <- seq(1, 100)
##' incidence <- rep(NA, length(day))
##' true_history <- array(NA_real_, c(5, 1, 101))
##' true_history[, 1, 1] <- sir$state()
##' for (i in day) {
##'   state_start <- sir$state()
##'   sir$run(i / dt)
##'   state_end <- sir$state()
##'   true_history[, 1, i + 1] <- state_end
##'   # Reduction in S
##'   incidence[i] <- state_start[1, 1] - state_end[1, 1]
##' }
##'
##' # Convert this into our required format:
##' data_raw <- data.frame(day = day, incidence = incidence)
##' data <- particle_filter_data(data_raw, "day", 4)
##'
##' # A comparison function
##' compare <- function(state, observed, pars = NULL) {
##'   if (is.null(pars$exp_noise)) {
##'     exp_noise <- 1e6
##'   } else {
##'     exp_noise <- pars$exp_noise
##'   }
##'   incidence_modelled <- state[1,]
##'   incidence_observed <- observed$incidence
##'   lambda <- incidence_modelled +
##'     rexp(length(incidence_modelled), exp_noise)
##'   dpois(incidence_observed, lambda, log = TRUE)
##' }
##'
##' # Range and initial values for model parameters
##' pars <- mcstate::if2_parameters$new(
##'   list(mcstate::if2_parameter("beta", 0.15, min = 0, max = 1),
##'        mcstate::if2_parameter("gamma", 0.05, min = 0, max = 1)))
##'
##' # Set up of IF2 algorithm (the iterations and n_par_sets should be
##' # increased here for any real use)
##' control <- mcstate::if2_control(
##'   pars_sd = list("beta" = 0.02, "gamma" = 0.02),
##'   iterations = 10,
##'   n_par_sets = 40,
##'   cooling_target = 0.5,
##'   progress = interactive())
##'
##' # Create a particle filter object
##' filter <- mcstate::particle_filter$new(data, gen, 1L, compare)
##'
##' # Then run the IF2
##' res <- mcstate::if2(pars, filter, control)
##'
##' # Get log-likelihood estimates from running a particle filter at
##' # each final parameter estimate
##' ll_samples <- mcstate::if2_sample(res, 20)
if2 <- function(pars, filter, control) {
  assert_is(pars, "if2_parameters")
  assert_is(filter, "particle_filter")
  assert_is(control, "if2_control")

  inputs <- filter$inputs()

  name_order <- match(pars$names(), names(control$pars_sd))
  if (any(is.na(name_order))) {
    missing <- pars$names()[is.na(name_order)]
    stop(sprintf("'%s' must be in control$pars_sd",
                 str_collapse(missing)), call. = FALSE)
  }
  pars_sd <- unlist(control$pars_sd[name_order])

  data_split <- df_to_list_of_lists(inputs$data)

  steps <- attr(inputs$data, "steps")
  n_steps <- nrow(steps)

  n_par_sets <- control$n_par_sets
  iterations <- control$iterations
  cooling_rate <- 1 / iterations
  alpha_cool <- control$cooling_target^cooling_rate

  pars_matrix <- pars$walk_initialise(n_par_sets, pars_sd)
  n_pars <- nrow(pars_matrix)

  model <- inputs$model$new(pars = pars$model(pars_matrix),
                            step = steps[[1L]],
                            n_particles = NULL,
                            n_threads = inputs$n_threads,
                            seed = inputs$seed,
                            pars_multi = TRUE)
  if (!is.null(inputs$index)) {
    model$set_index(inputs$index(model$info())$run)
  }

  ## NOTE: the [, 1L]..[[1L]] here assumes that observation parameters
  ## always shared, and are not sampled (i.e., same over
  ## simulation). We can't enforce that though.
  pars_compare <- pars$model(pars_matrix[, 1L, drop = FALSE])[[1L]]

  ## We'll collect into these:
  log_likelihood <- numeric(iterations)
  result_pars <- array(NA_real_, c(n_pars, n_par_sets, iterations))

  p <- pmcmc_progress(iterations, control$progress)

  for (i in seq_len(iterations)) {
    p()
    model$update_state(pars = pars$model(pars_matrix), step = steps[[1L]])
    for (t in seq_len(n_steps)) {
      step_end <- steps[t, 2L]
      state <- model$run(step_end)

      log_weights <- inputs$compare(state, data_split[[t]], pars_compare)

      if (!is.null(log_weights)) {
        weights <- scale_log_weights(log_weights + pars$prior(pars_matrix))
        log_likelihood[i] <- log_likelihood[i] + weights$average
        if (log_likelihood[i] == -Inf) {
          break
        }

        kappa <- particle_resample(weights$weights)
        model$reorder(kappa)
        pars_matrix <- pars$walk(pars_matrix[, kappa], pars_sd)
        model$update_state(pars = pars$model(pars_matrix),
                           set_initial_state = FALSE)
      }
    }
    result_pars[, , i] <- pars_matrix

    pars_sd <- pars_sd * alpha_cool
    pars_matrix <- pars$walk(pars_matrix, pars_sd)
  }

  ## pars dimensions are: n_pars, n_par_sets, iterations
  ## ll dimensions are: iterations, n_par_sets
  rownames(result_pars) <- pars$names()
  result <- list(log_likelihood = log_likelihood, pars = result_pars)

  ret <- list(result = result,
              pars = pars,
              control = control,
              filter = inputs)
  class(ret) <- "if2_fit"
  ret
}



##' @rdname if2
##' @param obj An object of class `if2_fit`, returned by `mcstate::if2()`
##'
##' @param n_particles The number of particles to simulate, for each
##'   IF2 parameter set
##' @export
if2_sample <- function(obj, n_particles) {
  assert_is(obj, "if2_fit")

  inputs <- obj$filter
  inputs$n_particles <- n_particles

  filter <- particle_filter_from_inputs(inputs)
  n_par_sets <- obj$control$n_par_sets

  ll <- numeric(n_par_sets)
  pars <- obj$pars$model(array_drop(
    obj$result$pars[, , obj$control$iterations, drop = FALSE], 3))

  p <- pmcmc_progress(n_par_sets, obj$control$progress)
  for (i in seq_len(n_par_sets)) {
    p()
    ll[i] <- filter$run(pars = pars[[i]])
  }

  ll
}
