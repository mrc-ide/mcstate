##' Run an IF2 inference.
##'
##' See:
##' Ionides EL, Nguyen D, Atchadé Y, Stoev S, King AA (2015).
##' "Inference for Dynamic and Latent Variable Models via Iterated,
##' Perturbed Bayes Maps."
##' PNAS, 112(3), 719–724. https://doi.org/10.1073/pnas.1410597112.
##'
##' @title Run iterated filtering (IF2 algorithm)
##'
##' @description Create an IF2 object for running
##'   and interacting with an IF2 inference.
##'
##' @export
##' @importFrom R6 R6Class
##' @examples
##' # A basic SIR model used in the particle filter example
##' dat <- example_sir()
##'
##' # Range and initial values for model parameters
##' pars <- if2_parameters$new(
##'           list(if2_parameter("beta", 0.15, min = 0, max = 1),
##'                if2_parameter("gamma", 0.05, min = 0, max = 1)))
##'
##' # Set up of IF2 algorithm
##' iterations <- 100
##' cooling_target <- 0.5
##' n_par_sets <- 300
##' control <- if2_control(pars_sd = list("beta" = 0.02, "gamma" = 0.02),
##'                         iterations = iterations,
##'                         n_par_sets = n_par_sets,
##'                         cooling_target = cooling_target,
##'                         progress = FALSE)
##'
##' Set up, and run IF2
##' filter <- if2$new(pars, dat$data, dat$model, dat$compare, NULL,
##'                   dat$index, control)
##' filter$run()
##'
##' # Plot results
##' filter$plot("ll")
##' filter$plot("beta")
##' filter$plot("gamma")
##'
##' # Get log-likelihood estimates from running a particle filter at
##' # each final parameter estimate
##' n_particles <- 100
##' ll_samples <- filter$sample(n_particles)
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

    ##' @description Create the IF2 object
    ##'
    ##' @param pars Parameters as an [if2_parameters()] object
    ##'
    ##' @param data The data set to be used, as in the [particle_filter()]
    ##'
    ##' @param model A stochastic model to use.  Must be a
    ##' `dust_generator` object.
    ##'
    ##' @param compare A comparison function, as in the [particle_filter()]
    ##'
    ##' @param compare_pars A list of any parameters required by the
    ##' comparison function
    ##'
    ##' @param index An index function, as in the [particle_filter()]
    ##'
    ##' @param control An [if2_control()] object
    ##'
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

    ##' @description Run the IF2
    ##'
    ##' @param n_threads Number of threads to use when running the
    ##' simulation. Defaults to 1, and should not be set higher than the
    ##' number of cores available to the machine.
    ##'
    ##' @param seed Seed for the random number generator on initial
    ##' creation. Can be `NULL` (to initialise using R's random number
    ##' generator), a positive integer, or a raw vector - see [`dust::dust`]
    ##' and [`dust::dust_rng`] for more details. Note that the random number
    ##' stream is unrelated from R's random number generator, except for
    ##' initialisation with `seed = NULL`.
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
        model$set_index(private$index(model$info())$run)
      }

      log_likelihood <- rep(0, iterations)
      if_pars <- array(NA_real_, c(n_pars, n_par_sets, iterations))
      alpha_cool <- cooling_target ^ (1 / iterations)

      p <- pmcmc_progress(iterations, private$control$progress)

      for (m in seq_len(iterations)) {
        p()
        model$reset(pars = private$pars$model(pars_matrix), steps[[1L]])
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
      # pars dimensions are: n_pars, n_par_sets, iterations
      # ll dimensions are: iterations, n_par_sets
      private$ll <- log_likelihood
      private$if_pars <- if_pars
    },

    ##' @description Return the log-likelihood at each iteration
    log_likelihood = function() {
      if (is.null(private$ll)) {
        stop("IF2 must be run first")
      }
      private$ll
    },

    ##' @description Return the set of parameters at each iteration
    pars_series = function() {
      if (is.null(private$ll)) {
        stop("IF2 must be run first")
      }
      private$if_pars
    },

    ##' @description Run a particle filter at each point estimate at final
    ##' state to get an estimate of the model likelihood
    ##'
    ##' @param n_particles The number of particles to simulate for each
    ##' parameter set
    ##'
    ##' @param n_threads Number of threads to use when running the
    ##' simulation. Defaults to 1, and should not be set higher than the
    ##' number of cores available to the machine.
    ##'
    ##' @param seed Seed for the random number generator on initial
    ##' creation. Can be `NULL` (to initialise using R's random number
    ##' generator), a positive integer, or a raw vector - see [`dust::dust`]
    ##' and [`dust::dust_rng`] for more details. Note that the random number
    ##' stream is unrelated from R's random number generator, except for
    ##' initialisation with `seed = NULL`.
    sample = function(n_particles, n_threads = 1L, seed = NULL) {
      if (is.null(private$ll)) {
        stop("IF2 must be run first")
      }

      n_par_sets <- private$control$n_par_sets
      n_iterations <- private$control$iterations
      pf_ll <- array(NA_real_, n_par_sets)
      pars_array <- private$if_pars[, , n_iterations]
      rownames(pars_array) <- private$pars$names()

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
          pf$run(pars = c(as.list(pars_array[, par_set]),
                          private$compare_pars))
      }
      pf_ll
    }
  )
)
