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
    index = NULL,
    control = NULL,
    pars_sd = NULL,
    # Outputs
    ll = NULL,
    if_pars = NULL
  ),

  public = list(

    initialize = function(pars, data, model, compare, compare_pars,
                          index, control) {
      assert_is(pars, "if2_parameters")
      assert_is(data, "particle_filter_data")
      if (!is_dust_generator(model)) {
        stop("'model' must be a dust_generator")
      }
      if (!is.null(index) && !is.function(index)) {
        stop("'index' must be function if not NULL")
      }
      if (!is.function(compare)) {
        stop("'compare' must be a function")
      }
      assert_is(control, "if2_control")
      name_order <- match(pars$names(), names(control$pars_sd))
      if (any(is.na(name_order))) {
        missing <- pars$names()[is.na(name_order)]
        stop(sprintf("'%s' must be in control$pars_sd",
                     str_collapse(missing)), call. = FALSE)
      }
      pars_sd <- unlist(control$pars_sd[name_order])

      private$pars <- pars
      private$data <- data
      private$model <- model
      private$compare <- compare
      private$compare_pars <- compare_pars
      private$index <- index
      private$control <- control
      private$pars_sd <- pars_sd
    },

    run = function(n_threads = 1L, seed = NULL) {
      data_split <- df_to_list_of_lists(private$data)

      steps <- unname(as.matrix(private$data[c("step_start", "step_end")]))
      n_steps <- nrow(steps)

      # Unpack some items from control
      n_par_sets <- private$control$n_par_sets
      iterations <- private$control$iterations
      cooling_target <- private$control$cooling_target
      pars_sd <- private$pars_sd

      pars_matrix <- private$pars$walk_initialise(n_par_sets,
                                                  pars_sd)
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

      p <- pmcmc_progress(iterations, private$control$progress)

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
            pars_matrix <- private$pars$walk(pars_matrix[, kappa],
                                             pars_sd)
            model$set_pars(private$pars$model(pars_matrix))
          }
        }
        pars_sd <- pars_sd * alpha_cool
        pars_final <- pars_matrix
        pars_matrix <- private$pars$walk(pars_final, pars_sd)
        if_pars[, , m] <- pars_final
      }
      # outputs
      # pars: n_pars * n_par_sets * iterations
      # ll: iterations * n_par_sets
      private$ll <- log_likelihood
      private$if_pars <- if_pars
    },

    ##' @description Return the initial parameter values as a named numeric
    ##' vector
    log_likelihood = function() {
      if(is.null(private$ll)) {
        stop("IF2 must be run first")
      }
      private$ll
    },

    pars_series = function() {
      if(is.null(private$ll)) {
        stop("IF2 must be run first")
      }
      private$if_pars
    },

    plot = function(what = "ll") {
      if(is.null(private$ll)) {
        stop("IF2 must be run first")
      }
      if (what %in% private$pars$names()) {
        par_idx <- which(private$pars$names() == what)
        mean <- apply(private$if_pars[par_idx, , ], 2, mean)
        quantiles <- apply(private$if_pars[par_idx, , ], 2,
                           quantile, c(0.025, 0.975))
        matplot(seq_len(length(private$ll)),
                mean, type = "l", lwd = 1, col = "#000000",
                xlab = "IF iteration", ylab = what,
                ylim = range(quantiles))
        matlines(seq_len(length(private$ll)),
                t(quantiles), type = "l", lty = 2, lwd = 1, col = "#999999")
        legend("bottomright", lwd = 1,
               legend = c("Mean", "95% quantile"), bty = "n")
      } else {
        plot(private$ll,
        main = "LL profile",
        xlab = "IF iteration",
        ylab = "log-likelihood",
        type = "l")
      }
    },

    # Run a particle filter at each point estimate at final state to get
    # mean + standard error
    sample = function(n_particles, n_threads = 1L, seed = NULL) {
      if(is.null(private$ll)) {
        stop("IF2 must be run first")
      }

      n_par_sets <- private$control$n_par_sets
      n_iterations <- private$control$iterations
      pf_ll <- array(NA_real_, n_par_sets)

      p <- pmcmc_progress(n_par_sets, private$control$progress)
      for (par_set in seq_len(n_par_sets)) {
        p()
        pf <- particle_filter$new(private$data,
                                  private$model,
                                  n_particles,
                                  private$compare,
                                  private$index,
                                  n_threads = n_threads,
                                  seed = seed)
        pf_ll[par_set] <-
          pf$run(pars = list(private$if_pars[, par_set, n_iterations]))
      }
      pf_ll
    }
  )
)
