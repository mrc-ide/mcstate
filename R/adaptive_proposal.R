##' Control for adaptive proposals, used in [mcstate::pmcmc_control]
##' for deterministic models.
##'
##' Efficient exploration of the parameter space during an MCMC might
##' be difficult when the target distribution is of high
##' dimensionality, especially if the target probability distribution
##' present a high degree of correlation.  Adaptive schemes are used
##' to "learn" on the fly the correlation structure by updating the
##' proposal distribution by recalculating the empirical
##' variance-covariance matrix and rescale it at each adaptive step of
##' the MCMC.
##'
##' Our implementation of an adaptive MCMC algorithm is based on an
##' adaptation of the "BlkAdpMul" algorithm in Sherlock et
##' al. (ALGORITHM 6). The algorithm is based on a random-walk
##' Metropolis-Hasting algorithm where the proposal is a multi-variate
##' Normal distribution centered on the current point.
##'
##' Sherlock C, Fearnhead P, Roberts GO (2010) The Random Walk
##' Metropolis: Linking Theory and Practice Through a Case
##' Study. Statistical Science 25:172–190.
##'
##' @title Adaptive proposal control
##'
##' @param initial_scaling The initial scaling of the variance
##'   covariance matrix to be used to generate the multivariate normal
##'   proposal for the random-walk Metropolis-Hasting algorithm.
##'
##' @param scaling_increment The scaling increment which is added or
##'   substracted to the scaling factor of the variance-covariance
##'   after each adaptive step.
##'
##' @param acceptance_target The target for the fraction of proposals
##'   that should be accepted (optimally) for the adaptive part of the
##'   mixture model.
##'
##' @param initial_weight Initial weight of the variance-covariance
##'   matrix used to build the proposal of the random-walk. Higher
##'   values translate into higher confidence of the initial
##'   variance-covariance matrix and means that update from additional
##'   samples will be slower.
##'
##' @param adaptive_contribution The fractional contribution of the
##'   adaptive part of the proposal. The proposal is based on a
##'   mixture model, with the non adaptive part used for the proposal
##'   with a probability 1-adaptive_contribution.
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
      assert_is(pars, "pmcmc_parameters")
      self$pars <- pars
      self$control <- control

      self$scaling <- control$initial_scaling
      self$weight <- control$initial_weight

      vcv <- pars$vcv()
      # self$zero <- apply(vcv == 0, 1, all)

      self$mean <- pars$mean()
      self$autocorrelation <- vcv +
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
        self$scaling <- max(self$scaling -
          self$control$acceptance_target * self$control$scaling_increment,
          self$control$scaling_increment)
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


adaptive_proposal_nested <- R6::R6Class(
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
      assert_is(pars, "pmcmc_parameters_nested")
      self$pars <- pars
      self$control <- control

      browser()

      n_varied <- length(pars$names("varied"))
      n_populations <- length(pars$populations())

      self$scaling <- list(fixed = control$initial_scaling,
                           varied = rep(control$initial_scaling, n_populations))
      self$weight <- list(fixed = control$initial_weight,
                          varied = rep(control$initial_weight, n_populations))

      vcv <- list(fixed = pars$vcv("fixed"),
                  varied = pars$vcv("varied"))

      self$mean <- list(fixed = pars$mean("fixed"),
                        varied = pars$mean("varied"))

      self$autocorrelation <- list(
        fixed = initial_autocorrelation(vcv$fixed, self$weight$fixed,
                                        self$mean$fixed),
        varied = Map(initial_autocorrelation, vcv$varied, self$weight$varied,
                     self$mean$varied))
    },

    propose = function(theta, type) {
      ## For varied, if one is adaptive all are adaptive
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

    update = function(theta, type, accept) {
      if (!self$proposal_was_adaptive) {
        return(invisible())
      }
      if (accept) {
        self$scaling <- self$scaling +
          (1 - self$control$acceptance_target) * self$control$scaling_increment
      } else {
        self$scaling <- max(self$scaling -
          self$control$acceptance_target * self$control$scaling_increment,
          self$control$scaling_increment)
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


initial_autocorrelation <- function(vcv, weight, mean) {
  vcv + weight / (weight - 1) * qp(mean)
}
