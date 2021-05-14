##' @export
if2 <- R6::R6Class(
  "if2",
  cloneable = FALSE,

  private = list(
    # Inputs
    pars = NULL,
    data = NULL,
    model = NULL,
    compare = NULL,
    compare_pars = NULL,
    index = NULL
    # Outputs
    log_likelihood = NULL,
    if_pars = NULL
  ),

  public = list(

    initialize = function(pars, data, model, compare, compare_pars, index) {
      assert_is(pars, "if2_parameters")
      assert_is(data, "particle_filter_data")
      if (!is_dust_generator(generator)) {
        stop("'model' must be a dust_generator")
      }
      if (!is.null(index) && !is.function(index)) {
        stop("'index' must be function if not NULL")
      }
      if (!is.function(comapare)) {
        stop("'compare' must be a function")
      }

      private$pars <- pars
      private$data <- data
      private$model <- model
      private$compare <- compare
      private$compare_pars <- compare_pars
      private$index <- index
    },

    run = function(pars_sd, iterations, n_par_sets, cooling_target,
                   n_threads = 1L, seed = NULL, progress = TRUE)
      data_split <- df_to_list_of_lists(private$data)

      steps <- unname(as.matrix(private$data[c("step_start", "step_end")]))
      n_steps <- nrow(steps)

      pars_matrix <- private$pars$walk_initialise(n_par_sets, pars_sd)
      n_pars <- nrow(pars_matrix)

      model <- private$model$new(pars = private$pars$model(pars_matrix),
                                 step = steps[[1L]],
                                 n_particles = NULL, n_threads = n_threads,
                                 seed = seed, pars_multi = TRUE)
      if (!is.null(private$index)) {
        model$set_index(private$index(mod$info())$run)
      }

      log_likelihood <- rep(0, iterations)
      if_pars <- array(NA_real_, c(n_pars, n_par_sets, iterations))
      alpha_cool <- cooling_target^(1 / iterations)

      p <- pmcmc_progress(iterations, progress)

      for (m in seq_len(iterations)) {
        p()
        model$reset(pars = private$pars$model(pars_matrix), steps[[1L]])
        ## TODO: refactor internal loop into function
        for (t in seq_len(n_steps)) {
          step_end <- steps[t, 2L]
          state <- model$run(step_end)

          log_weights <- private$compare(state, data_split[[t]],
                                         private$compare_pars)
          log_weights <- log_weights + private$pars$prior(pars_matrix)

          if (!is.null(log_weights)) {
            weights <- scale_log_weights(log_weights)
            log_likelihood[m] <- log_likelihood[m] + weights$average
            if (log_likelihood[m] == -Inf) {
              break
            }

            kappa <- particle_resample(weights$weights)
            model$reorder(kappa)
            pars_matrix <- pars$walk(pars_matrix[, kappa], pars_sd)
            model$set_pars(pars$model(pars_matrix))
          }
        }
        pars_sd <- pars_sd * alpha_cool
        pars_final <- pars_matrix
        pars_matrix <- pars$walk(pars_final, pars_sd)
        if_pars[, , m] <- pars_final
      }
      # outputs
      # pars: n_pars * n_par_sets * iterations
      # ll: iterations * n_par_sets
      private$log_likelihood <- log_likelihood
      private$if_pars <- if_pars
    },

    ##' @description Return the initial parameter values as a named numeric
    ##' vector
    log_likelihood = function() {
      private$log_likelihood
    },

    pars_series = function() {
      private$if_pars
    },

    ##' @description Set up a parameter walk
    ##'
    ##' @param n_par_sets An integer number of parameter sets, which
    ##' defines the size of the population being peturbed.
    walk_initialise = function(n_par_sets, pars_sd) {
      n_par_sets <- assert_integer(n_par_sets)
      n_pars <- length(private$parameters)

      pars_mat <- matrix(self$initial(), n_pars, n_par_sets)
      rownames(pars_mat) <- names(self$names())
      self$walk(pars_mat, pars_sd)
    },

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

  
}

##' @export
plot.if2_result <- function(obj, ...) {
  plot(obj$log_likelihood,
       main = "LL profile",
       xlab = "IF iteration",
       ylab = "log-likelihood",
       type = "l")

  # Could vectorise this better
  par_idx <- 0
  for (par in obj$pars$names()) {
    par_idx <- par_idx + 1
    mean <- apply(obj$if_pars[par_idx, , ], 2, mean)
    quantiles <- apply(obj$if_pars[par_idx, , ], 2, quantile, c(0.025, 0.975))
    matplot(seq_len(length(obj$log_likelihood)),
            mean, type = "l", lwd = 1, col = "#000000",
            xlab = "IF iteration", ylab = par,
            ylim = range(quantiles))
    matlines(seq_len(length(obj$log_likelihood)),
            t(quantiles), type = "l", lty = 2, lwd = 1, col = "#999999")
    legend("bottomright", lwd = 1, legend = c("Mean", "95% quantile"), bty = "n")
  }
}

# TODO add likelihood sample
# Run a particle filter at each point estimate at final state to get
# mean + standard error
summary
