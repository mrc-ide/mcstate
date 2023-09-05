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
##' @param diminishing_adaptation Logical parameter, set to `TRUE` as default.
##'   Indicates whether or not the adaptation of the scaling parameter
##'   diminishes with each update. Users should switch this option to `FALSE`
##'   with caution, as this will not ensure ergodicity of the MCMC algorithm.
##'   Thus we would not recommend production of posterior samples with this
##'   option set to `FALSE`, but users may want to use it for exploratory
##'   purposes.
##'
##' @export
adaptive_proposal_control <- function(initial_scaling = 0.2,
                                      scaling_increment = 0.01,
                                      acceptance_target = 0.234,
                                      initial_weight = 1000,
                                      adaptive_contribution = 0.95,
                                      diminishing_adaptation = TRUE) {
  ret <- list(initial_scaling = initial_scaling,
              scaling_increment = scaling_increment,
              acceptance_target = acceptance_target,
              initial_weight = initial_weight,
              adaptive_contribution = adaptive_contribution,
              diminishing_adaptation = diminishing_adaptation)
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

      self$mean <- pars$mean()
      self$autocorrelation <-
        initial_autocorrelation(vcv, self$weight, self$mean)
    },

    propose = function(theta) {
      self$proposal_was_adaptive <-
        runif(1) < self$control$adaptive_contribution
      vcv <- adaptive_vcv(self$scaling, self$autocorrelation, self$weight,
                          self$mean, self$proposal_was_adaptive)
      self$pars$propose(theta, vcv = vcv)
    },

    update = function(theta, accept) {
      if (!self$proposal_was_adaptive) {
        return(invisible())
      }
      self$weight <- self$weight + 1
      self$scaling <- update_scaling(self$scaling, self$weight,
                                     self$control, accept)
      self$autocorrelation <- update_autocorrelation(
        theta, self$weight, self$autocorrelation)
      self$mean <- update_mean(theta, self$weight, self$mean)
    }
  ))


adaptive_proposal_nested <- R6::R6Class(
  "adaptive_proposal",

  private = list(
    index = NULL
  ),

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

      nms <- pars$names("both")
      private$index <- list(fixed = match(pars$names("fixed"), nms),
                            varied = match(pars$names("varied"), nms))

      n_varied <- length(private$index$varied)
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
      if (type == "fixed") {
        self$proposal_was_adaptive <-
          runif(1) < self$control$adaptive_contribution
        vcv <- adaptive_vcv(self$scaling[[type]],
                            self$autocorrelation[[type]],
                            self$weight[[type]],
                            self$mean[[type]],
                            self$proposal_was_adaptive)
      } else if (type == "varied") {
        n_populations <- length(self$pars$populations())
        self$proposal_was_adaptive <-
          runif(n_populations) < self$control$adaptive_contribution
        vcv <- Map(adaptive_vcv, self$scaling[[type]],
                   self$autocorrelation[[type]],
                   self$weight[[type]],
                   self$mean[[type]],
                   self$proposal_was_adaptive)
      }
      self$pars$propose(theta, type, vcv = vcv)
    },

    update = function(theta, type, accept) {
      if (all(!self$proposal_was_adaptive)) {
        return(invisible())
      }

      idx <- private$index[[type]]
      if (type == "fixed") {
        theta_type <- theta[idx, 1, drop = TRUE]
      } else {
        theta_type <- lapply(seq_len(ncol(theta)), function(i)
                             theta[idx, i, drop = TRUE])
      }

      ## Probably we can save this more simply? - minor change on creation
      if (type == "fixed") {
        self$weight[[type]] <- self$weight[[type]] + 1
        self$scaling[[type]] <- update_scaling(
          self$scaling[[type]], self$weight[[type]], self$control, accept)
        self$autocorrelation[[type]] <- update_autocorrelation(
          theta_type, self$weight[[type]], self$autocorrelation[[type]])
        self$mean[[type]] <- update_mean(
          theta_type, self$weight[[type]], self$mean[[type]])
      } else if (type == "varied") {
        self$weight[[type]][self$proposal_was_adaptive] <-
          self$weight[[type]][self$proposal_was_adaptive] + 1
        self$scaling[[type]] <- update_scaling(
          self$scaling[[type]], self$weight[[type]], self$control, accept,
          self$proposal_was_adaptive)
        self$autocorrelation[[type]][] <- Map(
          update_autocorrelation, theta_type, self$weight[[type]],
          self$autocorrelation[[type]], self$proposal_was_adaptive)
        self$mean[[type]][] <- Map(
          update_mean, theta_type, self$weight[[type]], self$mean[[type]],
          self$proposal_was_adaptive)
      }
      

      ## This is where we're at now - there's an issue here with how
      ## we want to do the update - we can do it neatly enough with a
      ## Map below, but then we do need to if/else this over type
      ## which is ugly. Alternatively we might be able to nicely store
      ## things as vectors and matrices and recycle more natively but
      ## that affects the vcv step a little.
    }
  ))



qp <- function(x) {
  outer(x, x)
}


initial_autocorrelation <- function(vcv, weight, mean) {
  vcv + weight / (weight - 1) * qp(mean)
}


adaptive_vcv <- function(scaling, autocorrelation, weight, mean,
                         proposal_was_adaptive = TRUE) {
  if (proposal_was_adaptive) {
    vcv <- scaling * (autocorrelation - weight / (weight - 1) * qp(mean))
  } else {
    vcv <- NULL
  }
  vcv
}



update_scaling <- function(scaling, weight, control, accept,
                           proposal_was_adaptive = TRUE) {
  reject <- !accept
  acceptance_target <- control$acceptance_target
  scaling_increment <- control$scaling_increment
  if (control$diminishing_adaptation) {
    scale_inc <- scaling_increment / sqrt(weight - control$initial_weight)
  } else {
    scale_inc <- rep(scaling_increment, length(weight))
  }
  
  accept_update <- proposal_was_adaptive & accept
  if (any(accept_update)) {
    scaling[accept_update] <-
      (sqrt(scaling[accept_update]) +
         (1 - acceptance_target) * scale_inc[accept_update]) ^ 2
  }
  reject_update <- proposal_was_adaptive & reject
  if (any(reject_update)) {
    scaling[reject_update] <-
      pmax(sqrt(scaling[reject_update]) -
             acceptance_target * scale_inc[reject_update],
           scaling_increment) ^ 2
  }
  scaling
}


update_autocorrelation <- function(theta, weight, autocorrelation,
                                   proposal_was_adaptive = TRUE) {
  if (proposal_was_adaptive) {
    autocorrelation <-
      (1 - 1 / (weight - 1)) * autocorrelation + 1 / (weight - 1) * qp(theta)
  }
  autocorrelation
}


update_mean <- function(theta, weight, mean, proposal_was_adaptive = TRUE) {
  if (proposal_was_adaptive) {
    mean <- (1 - 1 / weight) * mean + 1 / weight * theta
  }
  mean
}
