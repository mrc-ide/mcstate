##' Control for adaptive proposals, used in [mcstate::pmcmc_control]
##' for deterministic models.
##'
##' @title Adaptive proposal control
##'
##' @param initial_scaling The initial scaling
##'
##' @param scaling_increment The scaling increment
##'
##' @param acceptance_target The fraction of proposals that should be
##'   accepted (optimally).
##'
##' @param initial_weight Initial weight
##'
##' @param adaptive_contribution The fractional contribution of the
##'   adaptive part of the proposal
##'
##' @export
adaptive_proposal_control <- function(initial_scaling = 0.2,
                                      scaling_increment = 0.01,
                                      acceptance_target = 0.234,
                                      initial_weight = 1000,
                                      adaptive_contribution = 0.95) {
  ret <- list(initial_scaling = initial_scaling,
              scaling_increment = scaling_increment,
              acceptance_target = acceptance_target,
              initial_weight = initial_weight,
              adaptive_contribution = adaptive_contribution)
  class(ret) <- "adaptive_proposal_control"
  ret
}


adaptive_proposal <- R6::R6Class(
  "adaptive_proposal",

  public = list(
    pars = NULL,
    control = NULL,

    mean = NULL,
    autocorrelation = NULL,
    weight = NULL,
    scaling = NULL,
    proposal_was_adaptive = NULL,

    initialize = function(pars, control) {
      self$pars <- pars
      self$control <- control

      self$scaling <- control$initial_scaling
      self$weight <- control$initial_weight

      self$mean <- pars$mean()
      self$autocorrelation <- pars$vcv() +
        self$weight / (self$weight - 1) * qp(self$mean)
    },

    propose = function(theta) {
      self$proposal_was_adaptive <-
        runif(1) < self$control$adaptive_contribution
      if (self$proposal_was_adaptive) {
        vcv <- self$scaling * (
          self$autocorrelation -
          self$weight / (self$weight - 1) * qp(self$mean))
        self$pars$propose(theta, vcv = vcv)
      } else {
        self$pars$propose(theta)
      }
    },

    update = function(theta, accept) {
      if (!self$proposal_was_adaptive) {
        return(invisible())
      }
      if (accept) {
        self$scaling <- self$scaling +
          (1 - self$control$acceptance_target) * self$control$scaling_increment
      } else {
        self$scaling <- self$scaling -
          self$control$acceptance_target * self$control$scaling_increment
      }

      ## Update of the autocorrelation matrix and mean of past samples
      self$weight <- self$weight + 1
      self$autocorrelation <-
        (1 - 1 / (self$weight - 1)) * self$autocorrelation +
        1 / (self$weight - 1) * qp(theta)
      self$mean <-
        (1 - 1 / self$weight) * self$mean +
        1 / self$weight * theta
    }
  ))


qp <- function(x) {
  outer(x, x)
}
