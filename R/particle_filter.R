##' @title Particle filter
##'
##' @description Create a \code{particle_filter} object for running
##'   and interacting with a particle filter.  A higher-level
##'   interface will be implemented later.
##'
##' @export
##' @examples
##' # A basic SIR model included in the package:
##' path <- system.file("example/sir/odin_sir.R", package = "mcstate")
##' gen <- odin::odin_(path)
##' sir <- gen()
##'
##' # Initial conditions for both the model and particle filter
##' y0 <- sir$initial(0)
##'
##' # Some data that we will fit to:
##' y <- sir$run(seq(0, 400, by = 4), y0)
##' data_raw <- as.data.frame(y)[c("day", "incidence")]
##'
##' # Convert this into our required format:
##' data <- mcstate::particle_filter_data(data_raw[-1, ], "day", 4)
##'
##' # A comparison function
##' compare <- function(state, output, observed, pars = NULL) {
##'   incid_modelled <- output[1, ]
##'   incid_observed <- observed$incidence
##'   lambda <- incid_modelled +
##'     rexp(n = length(incid_modelled), rate = 1e6)
##'   dpois(x = incid_observed, lambda = lambda, log = TRUE)
##' }
##'
##' # Construct the particle_filter object:
##' p <- mcstate::particle_filter$new(data, gen, compare)
##' p$run(y0, 100, TRUE)
##'
##' # Our simulated trajectories, with the "real" data superimposed
##' matplot(data_raw$day, t(p$history[1, , ]), type = "l",
##'         xlab = "Time", ylab = "State",
##'         col = "#ff000022", lty = 1, ylim = range(p$history))
##' matlines(data_raw$day, t(p$history[2, , ]), col = "#ffff0022", lty = 1)
##' matlines(data_raw$day, t(p$history[3, , ]), col = "#0000ff22", lty = 1)
##' matpoints(y[, "day"], y[, 2:4], pch = 19, col = c("red", "yellow", "blue"))
particle_filter <- R6::R6Class(
  "particle_filter",
  cloneable = FALSE,

  private = list(
    data = NULL,
    steps = NULL,
    n_steps = NULL,
    compare = NULL,
    last_model = NULL
  ),

  public = list(
    ##' @field model The dust model generator being simulated (cannot be
    ##' re-bound)
    model = NULL,

    ##' @field state The final state of the last run of the particle filter
    state = NULL,

    ##' @field history The history of the last run of the particle filter
    ##' (if enabled with \code{save_history = TRUE}, otherwise NULL
    history = NULL,

    ##' Create the particle filter
    ##'
    ##' @param data The data set to be used for the particle filter.
    ##' Must be a \code{\link{data.frame}} with at least columns
    ##' \code{step_start} and \code{step_end}.  Additional columns are
    ##' used for comparison with the simulation.
    ##'
    ##' @param model A stochastic model to use.  Must be an
    ##' \code{odin_generator} object (i.e., something that can be used to
    ##' create an \code{odin_model}).
    ##'
    ##' @param compare A comparison function.  Must take arguments
    ##' \code{state}, \code{output}, \code{data} and \code{pars} as arguments
    ##' (though the arguments may have different names).
    ##' \code{state} is the simulated model state (a matrix with as
    ##' many rows as there are state variables and as many columns as
    ##' there are particles.  \code{output} is the output variables, if
    ##' the model produces them (\code{NULL} otherwise) and \code{data}
    ##' is a \code{list} of observed data corresponding to the current
    ##' time's row in the \code{data} object provided here in the
    ##' constructor.  \code{pars} is any additional parameters passed
    ##' through to the comparison function (via the \code{pars_compare}
    ##' argument to \code{$run}).
    initialize = function(data, model, compare) {
      assert_is(attr(model, which="name", exact=TRUE), "dust_generator")
      self$model <- model
      private$data <- particle_filter_validate_data(data)
      private$steps <- cbind(vnapply(private$data, "[[", "step_start"),
                             vnapply(private$data, "[[", "step_end"))
      private$n_steps <- length(private$data)
      private$compare <- compare
      lockBinding("model", self)
    },

    ##' We probably need some special treatment for the initial case
    ##' but it's not clear that it belongs here, rather than in some
    ##' function above this, as state is just provided here as a vector
    ##'
    ##' Run the particle filter
    ##'
    ##' @param state The initial state. Can either be a vector (same
    ##' state for all particles) or a matrix with \code{n_particles}
    ##' columns
    ##'
    ##' @param n_particles The number of particles to simulate
    ##'
    ##' @param save_history Logical, indicating if the history of all
    ##' particles should be saved
    ##'
    ##' @param pars_model Optional parameters to use when creating the
    ##' model (\code{NULL} if the model should be default-constructed).
    ##'
    ##' @param pars_compare Optional parameters to use when applying
    ##' the \code{compare} function.  These parameters will be passed as
    ##' the 4th argument to \code{compare}.
    ##'
    ##' @param step_start Optional first step to start at.  If provided,
    ##' this must be within the range of the first epoch implied in your
    ##' \code{data} provided to the constructor (i.e., not less than the
    ##' first element of \code{step_start} and less than \code{step_end})
    ##'
    ##' @return A single numeric value representing the log-likelihood
    ##' (\code{-Inf} if the model is impossible)
    run = function(model_data, n_particles, save_history = FALSE,
                   pars_compare = NULL, step_start = NULL,
                   run_params) {

      compare <- private$compare
      steps <- particle_steps(private$steps, step_start)
      run_params <- validate_dust_params(run_params)

      model <- self$model$new(data = model_data, step = steps[1, 1], n_particles = n_particles,
                              n_threads = run_params["n_threads"],
                              n_generators = run_params["n_generators"],
                              seed = run_params["seed"])

      state <- model$state
      if (save_history) {
        history <- array(NA_real_, c(dim(state), private$n_steps + 1))
        history[, , 1] <- state
      } else {
        history <- NULL
      }

      log_likelihood <- 0
      for (t in seq_len(private$n_steps)) {
        model$run(steps[t, ])
        state <- model$state()
        log_weights <- compare(state, private$data[[t]], pars_compare)
        if (save_history) {
          history[, , t + 1L] <- state
        }

        if (!is.null(log_weights)) {
          weights <- scale_log_weights(log_weights)
          log_likelihood <- log_likelihood + weights$average
          if (weights$average == -Inf) {
            ## Everything is impossible, so stop here
            break
          }

          kappa <- particle_resample(weights$weights)
          model$reorder(kappa)
          if (save_history) {
            history <- history[, kappa, ]
          }
        }
      }

      self$history <- history
      self$state <- state
      private$last_model <- model

      log_likelihood
    },

    ##' Run the particle filter, having spread parameters out
    ##'
    ##' The interface here might change!
    ##'
    ##' @param state The initial state. Can either be a vector (same
    ##' state for all particles) or a matrix with \code{n_particles}
    ##' columns
    ##'
    ##' @param n_particles The number of particles to simulate
    ##'
    ##' @param pars A list of parameters
    ##'
    ##' @param index A parameter index
    ##'
    ##' @param save_history Logical, indicating if the history of all
    ##' particles should be saved
    run2 = function(n_particles, save_history = FALSE,
                   index, pars, run_params) {
      step_start <- pars_model <- pars_compare <- NULL
      if (length(index$step_start) > 0) {
        step_start <- pars[[index$step_start]]
      }

      if (length(index$pars_model) > 0) {
        model_data <- pars[index$model_data]
      }

      if (length(index$pars_compare) > 0) {
        pars_compare <- pars[index$pars_compare]
      }

      self$run(model_data, n_particles, save_history = FALSE,
               pars_compare, step_start = NULL, run_params)
    },

    ##' Create predicted trajectories, based on the final point of a
    ##' run with the particle filter
    ##'
    ##' @param t The steps to predict from, \emph{offset from the final
    ##' point}. As a result the first time-point of \code{t} must be 0.
    ##' The predictions will not however, include that point.
    ##'
    ##' @param append Logical, indicating if the predictions should be
    ##' appended onto the previous history of the simulation.
    ##'
    ##' @return A 3d array with dimensions representing (1) the state
    ##' vector, (2) the particle, (3) time
    predict = function(t, append = FALSE) {
      if (is.null(self$state)) {
        stop("Particle filter has not been run")
      }
      if (append && is.null(self$history)) {
        stop("Can't append without history")
      }
      if (t[[1]] != 0) {
        stop("Expected first 't' element to be zero")
      }
      step <- private$data[[private$n_steps]]$step_end + t

      res <- array(NA, dim = c(dim(self$state), length(step)))
      for (t in step) {
        private$last_model$run(step)
        res[,t] <- private$last_model$state()
      }
      if (append) {
        res <- dde_abind(self$history, res)
      }
      res
    }
  ))


particle_filter_validate_data <- function(data) {
  assert_is(data, "data.frame")
  msg <- setdiff(c("step_start", "step_end"), names(data))
  if (length(msg)) {
    stop("Expected columns missing from data: ",
         paste(squote(msg), collapse = ", "))
  }
  if (nrow(data) < 2) {
    stop("Expected at least two time windows")
  }
  ## TODO: step_start and step_end must be integer-like
  if (all(data$step_start[-1] != data$step_end[-nrow(data)])) {
    stop("Expected time windows to be adjacent")
  }

  ## Processing to make future use nicer:
  lapply(unname(split(data, seq_len(nrow(data)))), as.list)
}


##' @importFrom stats runif
particle_resample <- function(weights) {
  n <- length(weights)
  u <- runif(1, 0, 1 / n) + seq(0, by = 1 / n, length.out = n)
  cum_weights <- cumsum(weights / sum(weights))
  findInterval(u, cum_weights) + 1L
}


scale_log_weights <- function(log_weights) {
  max_log_weights <- max(log_weights)
  if (max_log_weights == -Inf) {
    ## if all log_weights at a time-step are -Inf, this should
    ## terminate the particle filter and output the marginal
    ## likelihood estimate as -Inf
    average <- -Inf
    weights <- rep(NaN, length(log_weights))
  } else {
    ## calculation of weights, there is some rescaling here to avoid
    ## issues where exp(log_weights) might give computationally zero
    ## values
    weights <- exp(log_weights - max_log_weights)
    average <- log(mean(weights)) + max_log_weights
  }
  list(weights = weights, average = average)
}


particle_initial_state <- function(state, n_particles) {
  if (is.matrix(state)) {
    if (ncol(state) != n_particles) {
      stop(sprintf("Expected '%d' columns for initial state", n_particles))
    }
  } else {
    state <- matrix(state, length(state), n_particles)
  }
  state
}


particle_steps <- function(steps, step_start) {
  if (!is.null(step_start)) {
    assert_integer(step_start)
    if (step_start < steps[1, 1, drop = TRUE]) {
      stop(sprintf(
        "'step_start' must be >= %d (the first value of data$step_start)",
        steps[1, 1, drop = TRUE]))
    }
    if (step_start >= steps[1, 2, drop = TRUE]) {
      stop(sprintf(
        "'step_start' must be < %d (the first value of data$step_end)",
        steps[1, 2, drop = TRUE]))
    }
    steps[1, 1] <- step_start
  }
  steps
}

validate_dust_params <- function(run_params) {
  if (run_params['n_threads'] == NULL ||
      !assert_integer(run_params['n_threads']) ||
      run_params['n_threads'] < 1) {
    run_params['n_threads'] <- 1
  }

  if (run_params['seed'] == NULL ||
      !assert_integer(run_params['seed'])) {
    run_params['seed'] <- 1
  }

  if (run_params['n_generators'] == NULL ||
      !assert_integer(run_params['n_generators']) ||
      run_params['n_generators'] < run_params['n_threads']) {
    run_params['n_generators'] <- run_params['n_threads']
  }
  run_params
}