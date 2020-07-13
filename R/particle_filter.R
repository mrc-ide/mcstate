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
##' history <- array(NA_real_, c(4, 1, 101))
##' history[, 1, 1] <- sir$state()
##' for (i in day) {
##'    state_start <- sir$state()
##'    sir$run(i / dt)
##'    state_end <- sir$state()
##'    history[, 1, i + 1] <- state_end
##'    # Reduction in S
##'    incidence[i] <- state_start[1, 1] - state_end[1, 1]
##' }
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
##' matpoints(data_raw$day, t(history[1:3, , -1]), pch = 19,
##'           col = c("red", "yellow", "blue"))
particle_filter <- R6::R6Class(
  "particle_filter",
  cloneable = FALSE,

  private = list(
    data = NULL,
    steps = NULL,
    index = NULL,
    initial = NULL,
    n_steps = NULL,
    compare = NULL,
    ## Updated when the model is run
    last_model = NULL,
    index_state = NULL
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
    ##' @param model A stochastic model to use.  Must be a
    ##' \code{dust_generator} object.
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
    ##'
    ##' @param index An index function. This is used to compute the
    ##' "interesting" indexes of your model. It must be a function of
    ##' one argument, which will be the result of calling the
    ##' \code{$info()} method on your model. It should return a list
    ##' with elements \code{run} (indices to return at the end of each
    ##' run, passed through to your compare function) and \code{state}
    ##' (indices to return if saving state). These indices can overlap
    ##' but do not have to. This argument is optional but using it will
    ##' likely speed up your simulation if you have more than a few
    ##' states as it will reduce the amount of memory copied back and
    ##' forth.
    ##'
    ##' @param initial A function to generate initial conditions. If
    ##' given, then this function must accept 3 arguments: \code{info}
    ##' (the result of calling \code{$info()} as for \code{index}),
    ##' \code{n_particles} (the number of particles that the particle
    ##' filter is using) and \code{pars} (parameters passed in in the
    ##' \code{$run} method via the \code{pars_initial} argument).  It
    ##' must return a list, which can have the elements \code{state}
    ##' (initial model state, passed to the particle filter - either a
    ##' vector or a matrix, and overriding the initial conditions
    ##' provided by your model) and \code{step} (the initial step,
    ##' overriding the first step of your data - this must occur within
    ##' your first epoch in your \code{data} provided to the
    ##' constructor, i.e., not less than the first element of
    ##' \code{step_start} and not more than \code{step_end})
    initialize = function(data, model, compare, index = NULL, initial = NULL) {
      if (!is_dust_generator(model)) {
        stop("'model' must be a dust_generator")
      }
      if (!is.null(index) && !is.function(index)) {
        stop("'index' must be function if not NULL")
      }
      if (!is.null(initial) && !is.function(initial)) {
        stop("'initial' must be function if not NULL")
      }

      self$model <- model
      private$data <- particle_filter_validate_data(data)
      private$steps <- cbind(vnapply(private$data, "[[", "step_start"),
                             vnapply(private$data, "[[", "step_end"))
      private$n_steps <- length(private$data)
      private$compare <- compare
      private$index <- index
      private$initial <- initial
      lockBinding("model", self)
    },

    ##' Run the particle filter
    ##'
    ##' @param pars_model The \code{data} object passed into dust,
    ##' which may contain parameters for your model (these might affect
    ##' running and/or initial conditions, depending on your model and
    ##' if you are using \code{initial} / \code{pars_initial} below).
    ##'
    ##' @param n_particles The number of particles to simulate
    ##'
    ##' @param save_history Logical, indicating if the history of all
    ##' particles should be saved
    ##'
    ##' @param pars_compare Optional parameters to use when applying
    ##' the \code{compare} function.  These parameters will be passed as
    ##' the 4th argument to \code{compare}.
    ##'
    ##' @param run_params List containing seed and n_threads for use with dust
    ##'
    ##' @param pars_initial Parameters passed through to the \code{initial}
    ##' function (if provided).
    ##'
    ##' @return A single numeric value representing the log-likelihood
    ##' (\code{-Inf} if the model is impossible)
    run = function(pars_model, n_particles, save_history = FALSE,
                   pars_compare = NULL, run_params = NULL,
                   pars_initial = NULL) {
      compare <- private$compare
      steps <- private$steps
      run_params <- validate_dust_params(run_params)

      if (is.null(private$last_model)) {
        model <- self$model$new(data = pars_model, step = steps[[1L]],
                                n_particles = n_particles,
                                n_threads = run_params[["n_threads"]],
                                seed = run_params[["seed"]])
      } else {
        model <- private$last_model
        model$reset(pars_model, steps[[1L]])
      }

      if (!is.null(private$initial)) {
        initial <- private$initial(model$info(), n_particles, pars_initial)
        if (is.list(initial)) {
          steps <- particle_steps(private$steps, initial$step)
          model$set_state(initial$state, initial$step)
        } else {
          model$set_state(initial)
        }
      } else if (!is.null(pars_initial)) {
        stop(paste("'pars_initial' was provided but you do not have",
                   "an 'initial' function to pass it to"))
      }

      if (is.null(private$index)) {
        index_state <- NULL
      } else {
        index <- private$index(model$info())
        if (!is.null(index$run)) {
          model$set_index(index$run)
        }
        index_state <- index$state
      }

      ## Baseline
      ##
      ## TODO(#22): This needs dealing with in the vignette (documenting
      ## that we need a dummy step here most likely if the user
      ## changes the initial location *and* if they need the history -
      ##
      ## TODO(dust#47): It would be nicer to have a better way of
      ## controlling this, especially to deal with irregular data.
      prev_res <- model$run(steps[[1]])

      unique_particles <- rep(n_particles, private$n_steps + 1)
      if (save_history) {
        state <- model$state(index_state)
        history <- array(NA_real_, c(dim(state), private$n_steps + 1))
        history[, , 1] <- state
      } else {
        history <- NULL
      }

      log_likelihood <- 0
      for (t in seq_len(private$n_steps)) {
        res <- model$run(steps[t, 2])
        if (save_history) {
          history[, , t + 1L] <- model$state(index_state)
        }

        log_weights <- compare(res, prev_res, private$data[[t]],
                               pars_compare)
        if (is.null(log_weights)) {
          prev_res <- res
        } else {
          weights <- scale_log_weights(log_weights)
          log_likelihood <- log_likelihood + weights$average
          if (weights$average == -Inf) {
            ## Everything is impossible, so stop here
            break
          }

          kappa <- particle_resample(weights$weights)
          unique_particles[t + 1L] <- length(unique(kappa))
          model$reorder(kappa)
          ## Because we will use the values from "run" a second time
          ## in the compare function, we need to reorder theem as
          ## well.
          prev_res <- res[, kappa, drop = FALSE]
          if (save_history) {
            ## TODO(#28): possible churn here
            history <- history[, kappa, , drop = FALSE]
          }
        }
      }

      self$history <- history
      self$unique_particles <- unique_particles
      self$state <- model$state(index_state)
      private$last_model <- model
      private$index_state <- index_state

      log_likelihood
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
    ##' @param run_params List containing seed and n_threads
    ##' for use with dust
    ##'
    run2 = function(n_particles, save_history = FALSE,
                    index, pars, run_params = NULL) {
      pars_model <- pars_compare <- pars_initial <- NULL
      if (length(index$pars_model) > 0) {
        pars_model <- pars[index$pars_model]
      }

      if (length(index$pars_compare) > 0) {
        pars_compare <- pars[index$pars_compare]
      }

      if (length(index$pars_initial) > 0) {
        pars_initial <- pars[index$pars_initial]
      }

      ## TODO(#27): run_params comes out and moves into the constructor
      self$run(pars_model, n_particles, save_history,
               pars_compare, run_params = run_params,
               pars_initial = pars_initial)
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
      model <- private$last_model

      step_end <- private$data[[private$n_steps]]$step_end
      if (model$step() != step_end) {
        ## TODO(#24): We need to be able to predict multiple times and that
        ## does not work as we lose the state here. This needs support
        ## in dust?
        ##
        ## If not we can reset the state and model time and state (that
        ## needs a little dust support too)
        stop("Can't yet run predict multiple times!")
      }
      step <- step_end + t

      index_state <- private$index_state
      res <- array(NA_real_, dim = c(dim(self$state), length(step) - 1))
      forecast_step <- 0L
      for (t in step[-1L]) {
        model$run(t)
        forecast_step <- forecast_step + 1L
        res[, , forecast_step] <- model$state(index_state)
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
  log_weights[is.nan(log_weights)] <- -Inf
  max_log_weights <- max(log_weights)
  if (!is.finite(max_log_weights)) {
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
    if (min(step_start) < steps[1, 1, drop = TRUE]) {
      stop(sprintf(
        "'step_start' must be >= %d (the first value of data$step_start)",
        steps[1, 1, drop = TRUE]))
    }
    if (max(step_start) > steps[1, 2, drop = TRUE]) {
      stop(sprintf(
        "'step_start' must be <= %d (the first value of data$step_end)",
        steps[1, 2, drop = TRUE]))
    }
    steps[1, 1] <- max(step_start)
  }
  steps
}


validate_dust_params <- function(run_params) {
  if (is.null(run_params)) {
    run_params <- list(n_threads = 1L, seed = 1L)
  } else {
    run_params[["n_threads"]] <-
      validate_dust_params_size(run_params[["n_threads"]])
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
