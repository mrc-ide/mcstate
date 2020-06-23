##' @title Particle filter
##'
##' @description Create a \code{particle_filter} object for running
##'   and interacting with a particle filter.  A higher-level
##'   interface will be implemented later.
##'
##' @export
##' @examples
##' # A basic SIR model included in the package:
##' path <- system.file("example/sir/dust_sir.cpp", package = "mcstate")
##' gen <- dust::dust(path)
##'
##' # Some data that we will fit to, using 1 particle:
##' sir <- gen$new(data = NULL, step = 0, n_particles = 1)
##' dt <- 1/4
##' day <- seq(1, 100)
##' incidence <- rep(NA, length(day))
##' history <- array(NA_real_, c(3, 1, 101))
##' history[, 1, 1] <- sir$state()
##' for (i in day) {
##'    state_start <- sir$state()
##'    sir$run(i / dt)
##'    state_end <- sir$state()
##'    history[, 1, i] <- state_end
##'    # Reduction in S
##'    incidence[i] <- state_start[1,1] - state_end[1,1]
##'  }
##'
##' # Convert this into our required format:
##' data_raw <- data.frame(day = day, incidence = incidence)
##' data <- particle_filter_data(data_raw, "day", 4)
##'
##' # A comparison function
##' compare <- function(state, prev_state, observed, pars = NULL) {
##'   if (is.null(pars)) {
##'     pars <- list(exp_noise = 1e6)
##'   }
##'   incidence_modelled <- prev_state[1,] - state[1,]
##'   incidence_observed <- observed$incidence
##'   lambda <- incidence_modelled +
##'     rexp(n = length(incidence_modelled), rate = pars$exp_noise)
##'   dpois(x = incidence_observed, lambda = lambda, log = TRUE)
##' }
##'
##' # Construct the particle_filter object:
##' p <- particle_filter$new(data, gen, compare)
##' p$run(NULL, 100, TRUE)
##'
##' # Our simulated trajectories, with the "real" data superimposed
##' matplot(data_raw$day, t(p$history[1, , -1]), type = "l",
##'          xlab = "Time", ylab = "State",
##'          col = "#ff000022", lty = 1, ylim = range(p$history))
##' matlines(data_raw$day, t(p$history[2, , -1]), col = "#ffff0022", lty = 1)
##' matlines(data_raw$day, t(p$history[3, , -1]), col = "#0000ff22", lty = 1)
##' matpoints(data_raw$day, t(history[, , -1]), pch = 19,
##'           col = c("red", "yellow", "blue"))
particle_filter <- R6::R6Class(
  "particle_filter",
  cloneable = FALSE,

  private = list(
    data = NULL,
    steps = NULL,
    particle_steps = NULL,
    current_step = NULL,
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

    log_likelihood = NULL,

    ##' @field history The history of the last run of the particle filter
    ##' (if enabled with \code{save_history = TRUE}, otherwise NULL
    history = NULL,

    ##' @field unique_particles The number of unique particles sampled
    ##' at each step that has been run
    unique_particles = NULL,

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
      if (!is_dust_generator(model)) {
        stop("'model' must be a dust_generator")
      }

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
    ##' @param model_data The data object passed into dust, which may contain
    ##' parameters and/or initial conditions
    ##'
    ##' @param n_particles The number of particles to simulate
    ##'
    ##' @param n_steps The number of steps to run. \code{NULL} runs to
    ##' the end of the data. (default = \code{NULL})
    ##'
    ##' @param save_history Logical, indicating if the history of all
    ##' particles should be saved
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
    ##' @param run_params List containing seed, n_threads and n_generators
    ##' for use with dust
    ##'
    ##' @return A single numeric value representing the log-likelihood
    ##' (\code{-Inf} if the model is impossible)
    run = function(model_data, n_particles, n_steps = NULL,
                   save_history = FALSE,
                   pars_compare = NULL, step_start = NULL,
                   run_params = NULL) {

      compare <- private$compare
      private$particle_steps <- particle_steps(private$steps, step_start)
      private$current_step <- 0
      run_params <- validate_dust_params(run_params)
      if (is.null(n_steps)) {
        n_steps <- private$n_steps
      } else if (n_steps > private$n_steps) {
        stop("Cannot run past end of data")
      }

      model <- self$model$new(data = model_data, step = steps[1, 1],
                              n_particles = n_particles,
                              n_threads = run_params[["n_threads"]],
                              n_generators = run_params[["n_generators"]],
                              seed = run_params[["seed"]])
      private$last_model <- model

      state <- model$state()
      self$unique_particles <- rep(n_particles, private$n_steps + 1)
      if (save_history) {
        history <- array(NA_real_, c(dim(state), private$n_steps + 1))
        history[, , 1] <- state
      } else {
        history <- NULL
      }
      self$history <- history

      log_likelihood <- 0
      for (t in seq_len(n_steps)) {
        state <- self$continue(state, pars_compare)
      }

      self$state <- state

      self$log_likelihood
    },

    ##' Run the particle filter, modifying parameters
    ##'
    ##' @param n_particles The number of particles to simulate
    ##'
    ##' @param save_history Logical, indicating if the history of all
    ##' particles should be saved
    ##'
    ##' @param index A parameter index
    ##'
    ##' @param pars A list of parameters
    ##'
    ##' @param run_params List containing seed, n_threads and n_generators
    ##' for use with dust
    ##'
    run2 = function(n_particles, n_steps = NULL,
                    save_history = FALSE,
                    index, pars, run_params = NULL) {
      step_start <- pars_model <- pars_compare <- NULL
      if (length(index$step_start) > 0) {
        step_start <- pars[[index$step_start]]
      }

      if (length(index$model_data) > 0) {
        model_data <- pars[index$model_data]
      }

      if (length(index$pars_compare) > 0) {
        pars_compare <- pars[index$pars_compare]
      }

      self$run(model_data, n_particles, n_steps, save_history,
               pars_compare, step_start = step_start, run_params = run_params)
    },

    # TODO - might be best if this didn't take any arguments
    # and data was just stored in the object
    continue = function(prev_state, pars_compare) {
      t = private$particle_steps[private$current_step, 2]
      private$current_step <- private$current_step + 1
      private$last_model$run(t)

      state <- model$state()
      log_weights <- compare(state, prev_state, private$data[[t]],
                              pars_compare)
      if (save_history) {
        private$history[, , t + 1L] <- state
      }

      if (!is.null(log_weights)) {
        weights <- scale_log_weights(log_weights)
        self$log_likelihood <- self$log_likelihood + weights$average
        if (weights$average == -Inf) {
          ## Everything is impossible, so stop here
          break
        }

        kappa <- particle_resample(weights$weights)
        self$unique_particles[t + 1L] <- length(unique(kappa))
        model$reorder(kappa)
        state <- state[, kappa]
        if (save_history) {
          private$history <- private$history[, kappa, ]
        }
        state
      }
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

      res <- array(NA_real_, dim = c(dim(self$state), length(step) - 1))
      forecast_step <- 0
      for (t in step[-1]) {
        private$last_model$run(t)
        forecast_step <- forecast_step + 1
        res[, , forecast_step] <- private$last_model$state()
      }
      if (append) {
        res <- dde_abind(self$history, res)
      }
      self$state <- private$last_model$state()
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
  if (is.null(run_params)) {
    run_params <- list(n_threads = 1L, n_generators = 1L, seed = 1L)
  } else {
    run_params[["n_threads"]] <-
      validate_dust_params_size(run_params[["n_threads"]])
    run_params[["n_generators"]] <-
      validate_dust_params_size(run_params[["n_generators"]])
    run_params[["seed"]] <-
      validate_dust_params_size(run_params[["seed"]])
  }
  run_params
}


validate_dust_params_size <- function(x) {
  if (is.null(x) || x < 1) {
    1L
  } else {
    as.integer(x)
  }
}


is_dust_generator <- function(x) {
  inherits(x, "R6ClassGenerator") &&
    identical(attr(x, which = "name", exact = TRUE), "dust_generator")
}
