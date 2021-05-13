##' Run a SMC^2. This is experimental and subject to change. Use at
##' your own risk.
##'
##' @title Run SMC^2
##'
##' @param pars A [`smc2_parameters`] object containing
##'   information about parameters (ranges, priors, proposal kernel,
##'   translation functions for use with the particle filter).
##'
##' @param filter A [`particle_filter`] object
##'
##' @param control A [mcstate::smc2_control] object to control the
##'   behaviour of the algorithm
##'
##' @return A `smc2_result` object, with elements
##'
##' * `pars`: a matrix of sampled parameters (n_parameter_set long)
##' * `probabilities`: a matrix of probabilities (`log_prior`, `log_likelihood`,
##'    `log_posterior` and `weight`). The latter is the log posterior
##'    normalised over all samples
##' * `statistics`: interesting or useful statistics about your sample,
##'    including the `ess` (effective sample size, over time),
##'    `acceptance_rate` (where a regeneration step was done, the acceptance
##'    rate), `n_particles`, `n_parameter_sets` and `n_steps` (inputs to
##'    the simulation). The `effort` field is a rough calculation of the
##'    number of particle-filter runs that this run was worth.
##' @export
##' @example man-roxygen/example-smc2.R
smc2 <- function(pars, filter, control) {
  assert_is(pars, "smc2_parameters")
  assert_is(filter, "particle_filter")
  assert_is(control, "smc2_control")
  obj <- smc2_engine$new(pars, filter, control)
  obj$run()
  obj$results()
}


smc2_engine <- R6::R6Class(
  "smc2_engine",

  private = list(
    pars = NULL,
    filter = NULL,
    state = NULL,
    n_steps = NULL,

    ## Monitoring
    acceptance_rate = NULL,
    ess = NULL,
    step_current = 0L,

    control = NULL,
    tick = NULL
  ),

  public = list(
    initialize = function(pars, filter, control) {
      private$pars <- pars
      private$filter <- filter
      private$control <- control

      theta <- pars$sample(control$n_parameter_sets)
      log_prior <- pars$prior(theta)
      log_likelihood <- rep(0.0, control$n_parameter_sets)

      inputs <- filter$inputs()
      ## TODO: fix this in dust to make it easier to get a
      ## "reasonable" state. We might also advance the state with a
      ## long jump first?
      seed <- dust::dust_rng$new(inputs$seed)$state()
      filters <- vector("list", control$n_parameter_sets)
      pars_model <- pars$model(theta)
      for (i in seq_len(control$n_parameter_sets)) {
        seed <- dust::dust_rng_state_long_jump(seed)
        f <- particle_filter$new(
          inputs$data, inputs$model, inputs$n_particles, inputs$compare,
          inputs$index, inputs$initial, inputs$n_threads, seed)
        filters[[i]] <- f$run_begin(pars_model[[i]], control$save_trajectories)
      }

      private$state <- smc2_state(filters, theta, log_prior, log_likelihood)

      private$n_steps <- nrow(inputs$data)
      private$acceptance_rate <- rep(NA_real_, private$n_steps)
      private$ess <- rep(NA_real_, private$n_steps)

      private$tick <- pmcmc_progress(private$n_steps, control$progress)
    },

    step = function() {
      t <- private$step_current + 1L
      step_ll <- vnapply(private$state$filter, function(f) f$step(t, TRUE))
      private$state$log_likelihood <- private$state$log_likelihood + step_ll
      private$state$log_posterior <- private$state$log_posterior + step_ll
      private$state$log_weights <- private$state$log_weights + step_ll
      private$state$weights <-
        scale_log_weights(private$state$log_weights)$weights
      ess <- sum(private$state$weights)^2 / sum(private$state$weights^2)
      private$step_current <- t
      private$ess[private$step_current] <- ess
      ess_crit <- private$control$degeneracy_threshold *
        length(private$state$filter)
      if (ess < ess_crit) {
        kernel <- self$vcv() * private$control$covariance_scaling
        self$resample()
        self$update(self$propose(kernel))
      }
      private$tick()
    },

    run = function() {
      for (i in seq_len(private$n_steps)) {
        self$step()
      }
    },

    vcv = function() {
      cov.wt(private$state$pars, wt = private$state$weights)$cov
    },

    resample = function() {
      ## Resample from the parameter sets according to their weights
      kappa <- particle_resample(private$state$log_weights)

      private$state$pars <- private$state$pars[kappa, ]
      private$state$log_prior <- private$state$log_prior[kappa]
      private$state$log_likelihood <- private$state$log_likelihood[kappa]
      private$state$log_posterior <- private$state$log_posterior[kappa]

      ## TODO: At this point we need to be ready to update our partial
      ## runs; after resampling if they've changed. This requires a
      ## full update of the state, unfortunately. We have to grab
      ## everything and do this all at once, but once working we can
      ## do this somewhat more space-efficiently if this causes too
      ## many problems by updating any index behind the current.
      filter_state <- lapply(private$state$filter, function(x)
        list(state = x$model$state(),
             log_likelihood = x$log_likelihood,
             history = x$history))
      for (i in seq_along(private$state$filter)) {
        s <- filter_state[[i]]
        private$state$filter[[i]]$model$set_state(s$state)
        private$state$filter[[i]]$log_likelihood <- s$log_likelihood
        private$state$filter[[i]]$history <- s$history
      }

      ## This should not really need doing as it will be reset later
      ## on if this has been called?
      private$state$log_weights <- private$state$log_weights[kappa]
    },

    propose = function(kernel) {
      pars <- private$pars$propose(private$state$pars, kernel)
      pars_model <- private$pars$model(pars)

      ## Calculate log prior for the proposed parameter sets
      log_prior <- private$pars$prior(pars)
      log_likelihood <- rep(-Inf, length(log_prior))
      filter <- vector("list", length(private$state$filter))

      ## NOTE: This could be done sequentially or concurrently. However,
      ## as we'll move to getting the particle filter set up to run with
      ## combinations of parameter sets soonish this is written to support
      ## that way.
      for (i in seq_along(filter)) {
        if (is.finite(log_prior[[i]])) {
          filter[[i]] <- private$state$filter[[i]]$fork(pars_model[[i]])
          log_likelihood[[i]] <- filter[[i]]$log_likelihood
        }
      }

      smc2_state(filter, pars, log_prior, log_likelihood)
    },

    update = function(proposed) {
      pr <- exp(proposed$log_posterior - private$state$log_posterior)
      accept <- runif(length(pr)) < pr
      reject <- !accept

      if (any(accept)) {
        private$state$pars[accept, ] <- proposed$pars[accept, ]
        private$state$log_prior[accept] <- proposed$log_prior[accept]
        private$state$log_likelihood[accept] <- proposed$log_likelihood[accept]
        private$state$log_posterior[accept] <- proposed$log_posterior[accept]
        private$state$filter[accept] <- proposed$filter[accept]
        private$state$log_weights[] <- 0
      }

      private$acceptance_rate[private$step_current] <- mean(accept)
    },

    results = function() {
      weight <- scale_log_weights(private$state$log_posterior)$weights
      probabilities <- cbind(
        log_prior = private$state$log_prior,
        log_likelihood = private$state$log_likelihood,
        log_posterior = private$state$log_posterior,
        weight = normalise(weight))
      pars <- private$state$pars

      ## The total number of filter steps we have taken; each time we
      ## resample we have to catch up *and* run the current set. We'd
      ## normalise this by dividing by n_steps (total number of full
      ## runs) and multiplying by n_parameter_sets (number of
      ## replicates) to get the equivalent number of complete runs of
      ## the particle filter (and that is roughly equal to the number
      ## of pmcmc steps)
      acc <- private$acceptance_rate
      effort <- sum(ifelse(is.na(acc), 1L, seq_along(acc) + 1L)) /
        private$n_steps * nrow(pars)

      ## TODO: Still don't get the trajectories out here
      ret <- list(pars = pars,
                  probabilities = probabilities,
                  statistics = list(
                    ess = private$ess,
                    acceptance_rate = private$acceptance_rate,
                    n_particles = private$filter$n_particles,
                    n_parameter_sets = nrow(pars),
                    n_steps = private$n_steps,
                    effort = effort))
      class(ret) <- "smc2_result"
      ret
    }
  ))


smc2_state <- function(filter, pars, log_prior, log_likelihood) {
  list(filter = filter,
       pars = pars,
       log_prior = log_prior,
       log_likelihood = log_likelihood,
       log_posterior = log_prior + log_likelihood,
       log_weights = rep(0, length(filter)),
       weights = rep(1, length(filter)))
}


##' @export
predict.smc2_result <- function(object, n = NULL, ...) {
  w <- object$probabilities[, "weight"]
  i <- sample.int(length(w), n %||% length(w), replace = TRUE, prob = w)
  object$pars[i, , drop = FALSE]
}
