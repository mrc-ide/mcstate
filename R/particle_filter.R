particle_filter <- R6::R6Class(
  "particle_filter",
  public = list(
    data = NULL,
    steps = NULL,
    n_steps = NULL,
    compare = NULL,
    save_history = NULL,

    state = NULL,
    history = NULL,

    initialize = function(data, compare, save_history = TRUE) {
      self$data <- particle_filter_validate_data(data)
      self$steps <- cbind(vnapply(self$data, "[[", "step_start"),
                          vnapply(self$data, "[[", "step_end"))
      self$n_steps <- length(self$data)
      self$compare <- compare
      self$save_history <- save_history
    },

    ## We probably need some special treatment for the initial case
    ## but it's not clear that it belongs here, rather than in some
    ## function above this, as state is just provided here as a vector
    run = function(state, model, n_particles) {
      if (is.matrix(state)) {
        browser()
      } else {
        state <- matrix(state, length(state), n_particles)
      }
      log_likelihood <- 0
      save_history <- self$save_history
      if (save_history) {
        history <- array(NA_real_, c(dim(state), self$n_steps + 1))
        history[, , 1] <- state
      } else {
        history <- NULL
      }
      for (t in seq_len(self$n_steps)) {
        res <- model$run(self$steps[t, ], state, replicate = n_particles,
                         use_names = FALSE, return_minimal = TRUE)
        state <- drop_dim(res, 2)
        output <- drop_dim(attr(res, "output", exact = TRUE), 2)
        if (save_history) {
          history[, , t + 1L] <- state
        }

        log_weights <- self$compare(state, output, self$data[[t]])

        if (!is.null(log_weights)) {
          weights <- scale_log_weights(log_weights)
          log_likelihood <- log_likelihood + weights$average
          if (weights$average == -Inf) {
            ## Everything is impossible, so stop here
            break
          }

          kappa <- particle_resample(weights$weights)
          state <- state[, kappa, drop = FALSE]
          if (save_history) {
            history <- history[, kappa, ]
          }
        }
      }

      self$history <- history
      self$state <- state

      log_likelihood
    },

    predict = function(t) {
      ## TODO: option of combined?
      step <- self$data[[self$n_steps]]$step_end + t
      n_particles <- ncol(self$state)
      self$model$run(step, self$state, replicate = n_particles,
                     use_names = FALSE, return_minimal = TRUE)
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
  if (max_log_weights == -Inf){
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