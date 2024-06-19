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
##' adaptation of the "accelerated shaping" algorithm in Spencer (2021).
##' The algorithm is based on a random-walk Metropolis-Hasting algorithm where
##' the proposal is a multi-variate Normal distribution centered on the current
##' point.
##'
##' Spencer SEF (2021) Accelerating adaptation in the adaptive 
##' Metropolis–Hastings random walk algorithm. Australian & New Zealand Journal
##' of Statistics 63:468-484.
##'
##' @title Adaptive proposal control
##'
##' @param initial_vcv_weight Weight of the initial variance-covariance
##'   matrix used to build the proposal of the random-walk. Higher
##'   values translate into higher confidence of the initial
##'   variance-covariance matrix and means that update from additional
##'   samples will be slower.
##'
##' @param initial_scaling The initial scaling of the variance
##'   covariance matrix to be used to generate the multivariate normal
##'   proposal for the random-walk Metropolis-Hastings algorithm. To generate
##'   the proposal matrix, the weighted variance covariance matrix is
##'   multiplied by the scaling parameter squared times 2.38^2 / n_pars (where
##'   n_pars is the number of fitted parameters). Thus, in a Gaussian target
##'   parameter space, the optimal scaling will be around 1.
##'
##' @param initial_scaling_weight The initial weight used in the scaling update.
##'   The scaling weight will increase after the first `pre_diminish`
##'   iterations, and as the scaling weight increases the adaptation of the
##'   scaling diminishes. If `NULL` (the default) the value is
##'   5 / (acceptance_target * (1 - acceptance_target)).
##'
##' @param min_scaling The minimum scaling of the variance covariance
##'   matrix to be used to generate the multivariate normal proposal
##'   for the random-walk Metropolis-Hastings algorithm.
##'
##' @param scaling_increment The scaling increment which is added or
##'   subtracted to the scaling factor of the variance-covariance
##'   after each adaptive step. If `NULL` (the default) then an optimal
##'   value will be calculated.
##'
##' @param log_scaling_update Logical, whether or not changes to the
##'   scaling parameter are made on the log-scale.
##'
##' @param acceptance_target The target for the fraction of proposals
##'   that should be accepted (optimally) for the adaptive part of the
##'   mixture model.
##'
##' @param forget_rate The rate of forgetting early parameter sets from the
##'   empirical variance-covariance matrix in the MCMC chains. For example,
##'   `forget_rate = 0.2` (the default) means that once in every 5th iterations
##'   we remove the earliest parameter set included, so would remove the 1st
##'   parameter set on the 5th update, the 2nd on the 10th update, and so
##'   on. Setting `forget_rate = 0` means early parameter sets are never
##'   forgotten.
##'
##' @param forget_end The final iteration at which early parameter sets can
##'   be forgotten. Setting `forget_rate = Inf` (the default) means that the
##'   forgetting mechanism continues throughout the chains. Forgetting early
##'   parameter sets becomes less useful once the chains have settled into the
##'   posterior mode, so this parameter might be set as an estimate of how long
##'   that would take.
##'
##' @param adapt_end The final iteration at which we can adapt the multivariate
##'   normal proposal. Thereafter the empirical variance-covariance matrix, its
##'   scaling and its weight remain fixed. This allows the adaptation to be
##'   switched off at a certain point to help ensure convergence of the chain.
##'
##' @param pre_diminish The number of updates before adaptation of the scaling
##'   parameter starts to diminish. Setting `pre_diminish = 0` means there is
##'   diminishing adaptation of the scaling parameter from the offset, while
##'   `pre_diminish = Inf` would mean there is never diminishing adaptation.
##'   Diminishing adaptation should help the scaling parameter to converge
##'   better, but while the chains find the location and scale of the posterior
##'   mode it might be useful to explore with it switched off.
##'
##' @export
adaptive_proposal_control <- function(initial_vcv_weight = 1000,
                                      initial_scaling = 1,
                                      initial_scaling_weight = NULL,
                                      min_scaling = 0,
                                      scaling_increment = NULL,
                                      log_scaling_update = TRUE,
                                      acceptance_target = 0.234,
                                      forget_rate = 0.2,
                                      forget_end = Inf,
                                      adapt_end = Inf,
                                      pre_diminish = 0) {
  ret <- list(initial_vcv_weight = initial_vcv_weight,
              initial_scaling = initial_scaling,
              initial_scaling_weight = initial_scaling_weight,
              min_scaling = min_scaling,
              scaling_increment = scaling_increment,
              log_scaling_update = log_scaling_update,
              acceptance_target = acceptance_target,
              forget_rate = forget_rate,
              forget_end = forget_end,
              adapt_end = adapt_end,
              pre_diminish = pre_diminish)
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
    initial_vcv = NULL,
    vcv = NULL,
    weight = NULL,
    iteration = NULL,
    included = NULL,
    scaling = NULL,
    scaling_increment = NULL,
    scaling_weight = NULL,
    
    initialize = function(pars, control) {
      assert_is(pars, "pmcmc_parameters")
      self$pars <- pars
      self$control <- control
      
      self$iteration <- 0
      self$weight <- 0
      
      self$mean <- pars$mean()
      self$initial_vcv <- pars$vcv()
      n_pars <- length(self$mean)
      self$autocorrelation <-
        matrix(0, n_pars, n_pars, dimnames = dimnames(self$initial_vcv))
      self$vcv <- update_vcv(self$mean, self$autocorrelation, self$weight)
      
      self$included <- integer()
      
      self$scaling <- control$initial_scaling
      self$scaling_increment <- control$scaling_increment %||%
        calc_scaling_increment(n_pars, control$acceptance_target,
                               control$log_scaling_update)
      self$scaling_weight <- control$initial_scaling_weight %||%
        5 / (control$acceptance_target * (1 - control$acceptance_target))
    },

    propose = function(theta) {
      proposal_vcv <-
        calc_proposal_vcv(self$scaling, self$vcv, self$weight,
                          self$initial_vcv, self$control$initial_vcv_weight)
      self$pars$propose(theta, vcv = proposal_vcv)
    },

    update = function(theta, accept_prob, theta_history, i) {
      self$iteration <- self$iteration + 1
      if (self$iteration > self$control$adapt_end) {
        return(invisible())
      }
      
      if (self$iteration > self$control$pre_diminish) {
        self$scaling_weight <- self$scaling_weight + 1
      }
      
      is_replacement <- check_replacement(self$iteration, self$control)
      if (is_replacement) {
        theta_remove <- theta_history[[self$included[1]]]
        self$included <- c(self$included[-1L], self$iteration)
      } else {
        theta_remove <- NULL
        self$included <- c(self$included, self$iteration)
        self$weight <- self$weight + 1 
      }
      
      
      self$scaling <-
        update_scaling(self$scaling, self$scaling_weight, accept_prob,
                       self$scaling_increment, self$control$min_scaling,
                       self$control$acceptance_target,
                       self$control$log_scaling_update)
      self$autocorrelation <- update_autocorrelation(
        theta, self$weight, self$autocorrelation, theta_remove)
      self$mean <- update_mean(theta, self$weight, self$mean, theta_remove)
      self$vcv <- update_vcv(self$mean, self$autocorrelation, self$weight)
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
    initial_vcv = NULL,
    vcv = NULL,
    weight = NULL,
    iteration = NULL,
    included = NULL,
    scaling = NULL,
    scaling_increment = NULL,
    scaling_weight = NULL,
    
    initialize = function(pars, control) {
      assert_is(pars, "pmcmc_parameters_nested")
      self$pars <- pars
      
      if (length(control$initial_vcv_weight) == 1) {
        control$initial_vcv_weight <-
          list(fixed = control$initial_vcv_weight,
               varied = control$initial_vcv_weight)
      }
      self$control <- control

      nms <- pars$names("both")
      private$index <- list(fixed = match(pars$names("fixed"), nms),
                            varied = match(pars$names("varied"), nms))

      n_fixed <- length(private$index$fixed)
      n_varied <- length(private$index$varied)
      n_populations <- length(pars$populations())

      self$weight <- list(fixed = 0,
                          varied = 0)
      
      self$iteration <- list(fixed = 0,
                             varied = 0)
      self$included <- list(fixed = integer(),
                            varied = integer())
      
      self$mean <- list(fixed = pars$mean("fixed"),
                        varied = pars$mean("varied"))
      self$initial_vcv <- list(fixed = pars$vcv("fixed"),
                               varied = pars$vcv("varied"))
      self$autocorrelation <- list(
        fixed = matrix(0, n_fixed, n_fixed,
                       dimnames = dimnames(self$initial_vcv$fixed)),
        varied = lapply(self$initial_vcv$varied,
                        function(x) {matrix(0, n_varied, n_varied,
                                            dimnames = dimnames(x))}))
      self$vcv <- list(
        fixed = update_vcv(self$mean$fixed, self$autocorrelation$fixed,
                           self$weight$fixed),
        varied = Map(update_vcv, self$mean$varied, self$autocorrelation$varied,
                     self$weight$varied))
      
      self$scaling <- list(fixed = control$initial_scaling,
                           varied = rep(control$initial_scaling, n_populations))
      self$scaling_increment <- list(
        fixed = control$scaling_increment %||% 
          calc_scaling_increment(n_fixed, control$acceptance_target,
                                 control$log_scaling_update),
        varied = control$scaling_increment %||% 
          calc_scaling_increment(n_varied, control$acceptance_target,
                                 control$log_scaling_update)
      )
      
      scaling_weight <- control$initial_scaling_weight %||%
        5 / (control$acceptance_target * (1 - control$acceptance_target))
      self$scaling_weight <- list(fixed = scaling_weight,
                                  varied = scaling_weight)

    },

    propose = function(theta, type) {
      if (type == "fixed") {
        proposal_vcv <-
          calc_proposal_vcv(self$scaling[[type]], self$vcv[[type]],
                            self$weight[[type]], self$initial_vcv[[type]],
                            self$control$initial_vcv_weight[[type]])
      } else if (type == "varied") {
        proposal_vcv <-
          Map(calc_proposal_vcv, self$scaling[[type]], self$vcv[[type]],
              self$weight[[type]], self$initial_vcv[[type]],
              self$control$initial_vcv_weight[[type]])
      }
      self$pars$propose(theta, type, vcv = proposal_vcv)
    },

    update = function(theta, type, accept_prob, theta_history, i) {
      idx <- private$index[[type]]
      if (type == "fixed") {
        theta_type <- theta[idx, 1, drop = TRUE]
      } else {
        theta_type <- lapply(seq_len(ncol(theta)), function(j)
                             theta[idx, j, drop = TRUE])
      }
      
      self$iteration[[type]] <- self$iteration[[type]] + 1
      if (self$iteration[[type]] > self$control$adapt_end) {
        return(invisible())
      }
      
      if (self$iteration[[type]] > self$control$pre_diminish) {
        self$scaling_weight[[type]] <- self$scaling_weight[[type]] + 1
      }
      
      is_replacement <- 
        check_replacement(self$iteration[[type]], self$control)
      if (is_replacement) {
        theta_remove <- theta_history[[self$included[[type]][1]]]
        if (type == "fixed") {
          theta_remove <- theta_remove[idx, 1, drop = TRUE]
        } else {
          theta_remove <- lapply(seq_len(ncol(theta_remove)), function(j)
            theta_remove[idx, j, drop = TRUE])
        }
        self$included[[type]] <-
          c(self$included[[type]][-1L], i)
      } else {
        if (type == "fixed") {
          theta_remove <- NULL
        } else {
          theta_remove <- rep(list(NULL), length(theta_type))
        }
        self$included[[type]] <- c(self$included[[type]], i)
        self$weight[[type]] <- self$weight[[type]] + 1
      }
      
      if (type == "fixed") {
        self$scaling[[type]] <-
          update_scaling(self$scaling[[type]], self$scaling_weight[[type]],
                         accept_prob, self$scaling_increment[[type]],
                         self$control$min_scaling,
                         self$control$acceptance_target,
                         self$control$log_scaling_update)
        self$autocorrelation[[type]] <- update_autocorrelation(
          theta_type, self$weight[[type]], self$autocorrelation[[type]],
          theta_remove)
        self$mean[[type]] <- update_mean(
          theta_type, self$weight[[type]], self$mean[[type]], theta_remove)
        self$vcv[[type]] <- update_vcv(
          self$mean[[type]], self$autocorrelation[[type]], self$weight[[type]])
      } else if (type == "varied") {
        self$scaling[[type]] <-
          update_scaling(self$scaling[[type]], self$scaling_weight[[type]],
                         accept_prob, self$scaling_increment[[type]],
                         self$control$min_scaling,
                         self$control$acceptance_target,
                         self$control$log_scaling_update)
        self$autocorrelation[[type]][] <- Map(
          update_autocorrelation, theta_type, self$weight[[type]],
          self$autocorrelation[[type]], theta_remove)
        self$mean[[type]][] <- Map(
          update_mean, theta_type, self$weight[[type]], self$mean[[type]],
          theta_remove)
        self$vcv[[type]] <-
          Map(update_vcv, self$mean[[type]], self$autocorrelation[[type]],
              self$weight[[type]])
      }

    }
  ))



qp <- function(x) {
  outer(x, x)
}


calc_scaling_increment <- function(n_pars, acceptance_target,
                                   log_scaling_update) {
  if (log_scaling_update) {
    A <- -qnorm(acceptance_target / 2)

    scaling_increment <-
      (1 - 1 / n_pars) * (sqrt(2 * pi) * exp(A^2 / 2)) / (2 * A) +
      1 / (n_pars * acceptance_target * (1 - acceptance_target))
  } else {
    scaling_increment <- 1 / 100
  }

  scaling_increment
}


calc_proposal_vcv <- function(scaling, vcv, weight, initial_vcv,
                              initial_vcv_weight) {
  n_pars <- nrow(vcv)

  if (weight == 0) {
    weighted_vcv <- initial_vcv
  } else {
    weighted_vcv <-
      ((weight - 1) * vcv + (initial_vcv_weight + n_pars + 1) * initial_vcv) /
      (weight + initial_vcv_weight + n_pars + 1)
  }


  2.38^2 / n_pars * scaling^2 * weighted_vcv
}


check_replacement <- function(iteration, control) {
  is_forget_step <- floor(control$forget_rate * iteration) >
    floor(control$forget_rate * (iteration - 1))
  is_before_forget_end <- iteration <= control$forget_end

  is_forget_step & is_before_forget_end
}


update_scaling <- function(scaling, scaling_weight, accept_prob,
                           scaling_increment, min_scaling,
                           acceptance_target, log_scaling_update) {
  scaling_change <- scaling_increment * (accept_prob - acceptance_target) /
    sqrt(scaling_weight)

  if (log_scaling_update) {
    pmax(min_scaling, scaling * exp(scaling_change))
  } else {
    pmax(min_scaling, scaling + scaling_change)
  }

}


update_autocorrelation <- function(theta, weight, autocorrelation,
                                   theta_remove) {
  if (!is.null(theta_remove)) {
    if (weight > 2) {
      autocorrelation <-
        autocorrelation + 1 / (weight - 1) * (qp(theta) - qp(theta_remove))
    } else {
      autocorrelation <- autocorrelation + qp(theta) - qp(theta_remove)
    }
  } else {
    if (weight > 2) {
      autocorrelation <-
        (1 - 1 / (weight - 1)) * autocorrelation + 1 / (weight - 1) * qp(theta)
    } else {
      autocorrelation <- autocorrelation + qp(theta)
    }
  }

  autocorrelation
}


update_mean <- function(theta, weight, mean, theta_remove) {
  if (!is.null(theta_remove)) {
    mean <- mean + 1 / weight * (theta - theta_remove)
  } else {
    mean <- (1 - 1 / weight) * mean + 1 / weight * theta
  }
  mean
}


update_vcv <- function(mean, autocorrelation, weight) {
  if (weight > 1) {
    vcv <- autocorrelation - weight / (weight - 1) * qp(mean)
  } else {
    vcv <- 0 * autocorrelation
  }

  vcv
}
