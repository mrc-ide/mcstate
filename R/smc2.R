smc2_engine <- R6::R6Class(
  "smc2_engine",

  private = list(
    parameters = NULL,
    filter = NULL,
    state = NULL,
    n_steps = NULL,

    ## Monitoring
    acceptance_rate = NULL,
    ess = NULL,
    step_current = 0L,

    ## Configuration
    degeneracy_threshold = NULL,
    covariance_scaling = NULL
  ),

  public = list(
    initialize = function(parameters, filter, n_parameter_sets,
                          degeneracy_threshold = 0.5,
                          covariance_scaling = 0.5,
                          save_history = FALSE) {
      private$parameters <- parameters
      private$filter <- filter

      pars <- parameters$sample(n_parameter_sets)
      log_prior <- parameters$prior(pars)
      log_likelihood <- rep(0.0, n_parameter_sets)

      inputs <- filter$inputs()
      ## TODO: fix this in dust to make it easier to get a
      ## "reasonable" state. We might also advance the state with a
      ## long jump first?
      seed <- dust::dust_rng$new(inputs$seed)$state()
      filters <- vector("list", n_parameter_sets)
      pars_model <- parameters$model(pars)
      for (i in seq_len(n_parameter_sets)) {
        seed <- dust::dust_rng_state_long_jump(seed)
        f <- particle_filter$new(
          inputs$data, inputs$model, inputs$n_particles, inputs$compare,
          inputs$index, inputs$initial, inputs$n_threads, seed)
        ## TODO: Do we need to save the underlying filter objects? If
        ## we're going to do projection they might be needed, but
        ## that's going to be hard to swing potentially?
        filters[[i]] <- f$run_begin(pars_model[[i]], save_history)
      }

      private$state <- smc2_state(filters, pars, log_prior, log_likelihood)

      private$n_steps <- nrow(inputs$data)
      private$acceptance_rate <- rep(NA_real_, private$n_steps)
      private$ess <- rep(NA_real_, private$n_steps)

      private$degeneracy_threshold <- degeneracy_threshold
      private$covariance_scaling <- covariance_scaling
    },

    step = function() {
      step_ll <- vnapply(private$state$filter, function(f) f$step(TRUE))
      private$state$log_likelihood <- private$state$log_likelihood + step_ll
      private$state$log_posterior <- private$state$log_posterior + step_ll
      private$state$log_weights <- private$state$log_weights + step_ll
      private$state$weights <-
        scale_log_weights(private$state$log_weights)$weights
      ess <- sum(private$state$weights)^2 / sum(private$state$weights^2)
      private$step_current <- private$step_current + 1L
      private$ess[private$step_current] <- ess
      if (ess < private$degeneracy_threshold * length(private$state$filter)) {
        message("Degenerate - rejuvenating filters")
        kernel <- self$vcv() * private$covariance_scaling
        self$resample()
        self$update(self$propose(kernel))
      }
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
      pars <- private$parameters$propose(private$state$pars, kernel)
      pars_model <- private$parameters$model(pars)

      ## Calculate log prior for the proposed parameter sets
      log_prior <- private$parameters$prior(pars)
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
