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
##' @param initial_scaling The initial scaling of the variance
##'   covariance matrix to be used to generate the multivariate normal
##'   proposal for the random-walk Metropolis-Hasting algorithm.
##'
##' @param scaling_increment The scaling increment which is added or
##'   subtracted to the scaling factor of the variance-covariance
##'   after each adaptive step.
##'
##' @param acceptance_target The target for the fraction of proposals
##'   that should be accepted (optimally) for the adaptive part of the
##'   mixture model.
##'
##' @param initial_vcv_weight Initial weight of the variance-covariance
##'   matrix used to build the proposal of the random-walk. Higher
##'   values translate into higher confidence of the initial
##'   variance-covariance matrix and means that update from additional
##'   samples will be slower.
##'   
##' @param forget_rate The rate of forgetting early parameter sets from the 
##'   empirical variance-covariance matrix in the MCMC chains. For example,
##'   `forget_rate = 0.2` (the default) means that once in every 5th updates
##'   we remove the earliest parameter set included, so would remove the 1st
##'   parameter set on the 5th update, the 2nd on the 10th update, and so
##'   on. Setting `forget_rate = 0` means early parameter sets are never
##'   forgotten.
##'   
##' @param forget_end The final update number at which early parameter sets can
##'   be forgotten. Setting `forget_rate = Inf` (the default) means that the
##'   forgetting mechanism continues throughout the chains. Forgetting early
##'   parameter sets becomes less useful once the chains have settled into the
##'   posterior mode, so this parameter might be set as an estimate of how long
##'   that would take
##'   
##' @param pre_diminish The number of updates before adaptation of the scaling
##'   parameter starts to diminish. Setting `pre_diminish = 0` means there is
##'   diminishing adaptation of the scaling parameter from the offset, while
##'   `pre_diminish = Inf` would mean there is never diminishing adaptation.
##'   Diminishing adaptation should help the scaling parameter to converge
##'   better, but while the chains find the location and scale of the posterior
##'   mode it might be useful to explore with it switched off
##'
##' @export
adaptive_proposal_control <- function(initial_scaling = 0.2,
                                      scaling_increment = 0.01,
                                      acceptance_target = 0.234,
                                      initial_vcv_weight = 1000,
                                      forget_rate = 0.2,
                                      forget_end = Inf,
                                      pre_diminish = 0) {
  ret <- list(initial_scaling = initial_scaling,
              scaling_increment = scaling_increment,
              acceptance_target = acceptance_target,
              initial_vcv_weight = initial_vcv_weight,
              forget_rate = forget_rate,
              forget_end = forget_end,
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
    weight = NULL,
    scaling = NULL,
    iteration = NULL,
    included = NULL,
    
    initialize = function(pars, control) {
      assert_is(pars, "pmcmc_parameters")
      self$pars <- pars
      self$control <- control

      self$scaling <- control$initial_scaling
      
      vcv <- pars$vcv()

      self$iteration <- 0
      self$weight <- 0
      self$mean <- pars$mean()
      self$autocorrelation <- 0 * vcv
      self$included <- c()
    },

    propose = function(theta) {
      vcv <- adaptive_vcv(self$scaling, self$autocorrelation, self$weight,
                          self$mean, self$pars$vcv(),
                          self$control$initial_vcv_weight)
      self$pars$propose(theta, vcv = vcv)
    },

    update = function(theta, accept, theta_history, i) {
      self$iteration <- self$iteration + 1
      is_replacement <- check_replacement(self$iteration, self$control)
      if (is_replacement) {
        theta_remove <- theta_history[[self$included[1]]]
      } else {
        theta_remove <- NULL
        self$weight <- self$weight + 1 
      }
      
      self$scaling <- update_scaling(self$scaling, self$iteration,
                                     self$control, accept)
      self$autocorrelation <- update_autocorrelation(
        theta, self$weight, self$autocorrelation, theta_remove)
      self$mean <- update_mean(theta, self$weight, self$mean, theta_remove)
      self$included <- update_included(self$included, i, is_replacement)
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
    iteration = NULL,
    included = NULL,

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

      n_varied <- length(private$index$varied)
      n_populations <- length(pars$populations())

      self$scaling <- list(fixed = control$initial_scaling,
                           varied = rep(control$initial_scaling, n_populations))
      self$weight <- list(fixed = 0,
                          varied = 0)
      
      self$iteration <- list(fixed = 0,
                             varied = 0)
      self$included <- list(fixed = c(),
                            varied = c())

      vcv <- list(fixed = pars$vcv("fixed"),
                  varied = pars$vcv("varied"))
      self$mean <- list(fixed = pars$mean("fixed"),
                        varied = pars$mean("varied"))

      self$autocorrelation <- list(
        fixed = 0 * vcv$fixed,
        varied = Map(function (x) 0 * x, vcv$varied))
    },

    propose = function(theta, type) {
      if (type == "fixed") {
        vcv <- adaptive_vcv(self$scaling[[type]],
                            self$autocorrelation[[type]],
                            self$weight[[type]],
                            self$mean[[type]],
                            self$pars$vcv("fixed"),
                            self$control$initial_vcv_weight[[type]])
      } else if (type == "varied") {
        vcv <- Map(adaptive_vcv, self$scaling[[type]],
                   self$autocorrelation[[type]],
                   self$weight[[type]],
                   self$mean[[type]],
                   self$pars$vcv("varied"),
                   self$control$initial_vcv_weight[[type]])
      }
      self$pars$propose(theta, type, vcv = vcv)
    },

    update = function(theta, type, accept, theta_history, i) {
      
      idx <- private$index[[type]]
      if (type == "fixed") {
        theta_type <- theta[idx, 1, drop = TRUE]
      } else {
        theta_type <- lapply(seq_len(ncol(theta)), function(j)
                             theta[idx, j, drop = TRUE])
      }
      
      ## Probably we can save this more simply? - minor change on creation
      self$iteration[[type]] <- self$iteration[[type]] + 1
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
      } else {
        self$weight[[type]] <- self$weight[[type]] + 1
        if (type == "fixed") {
          theta_remove <- NULL
        } else {
          theta_remove <- rep(list(NULL), length(theta_type))
        }
      }
      
      
      if (type == "fixed") {
        self$scaling[[type]] <- update_scaling(
          self$scaling[[type]], self$iteration[[type]], self$control, accept)
        self$autocorrelation[[type]] <- update_autocorrelation(
          theta_type, self$weight[[type]], self$autocorrelation[[type]],
          theta_remove)
        self$mean[[type]] <- update_mean(
          theta_type, self$weight[[type]], self$mean[[type]], theta_remove)
        self$included[[type]] <- 
          update_included(self$included[[type]], i, is_replacement)
      } else if (type == "varied") {
        self$scaling[[type]] <- update_scaling(
          self$scaling[[type]], self$iteration[[type]], self$control, accept)
        self$autocorrelation[[type]][] <- Map(
          update_autocorrelation, theta_type, self$weight[[type]],
          self$autocorrelation[[type]], theta_remove)
        self$mean[[type]][] <- Map(
          update_mean, theta_type, self$weight[[type]], self$mean[[type]],
          theta_remove)
        self$included[[type]] <- 
          update_included(self$included[[type]], i, is_replacement)
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

check_replacement <- function(iteration, control) {
  is_forget_step <- floor(control$forget_rate * iteration) >
    floor(control$forget_rate * (iteration - 1))
  is_before_forget_end <- iteration <= control$forget_end
  
  is_forget_step & is_before_forget_end
}

adaptive_vcv <- function(scaling, autocorrelation, weight, mean, initial_vcv,
                         initial_vcv_weight) {
  if (weight > 1) {
    vcv <- autocorrelation - weight / (weight - 1) * qp(mean)
  } else {
    vcv <- 0 * autocorrelation
  }
  
  d <- length(mean)
  
  scaling * ((weight - 1) * vcv + (initial_vcv_weight + d + 1) * initial_vcv) /
    (weight + initial_vcv_weight + d + 1) 
    
}



update_scaling <- function(scaling, iteration, control, accept) {
  reject <- !accept
  acceptance_target <- control$acceptance_target
  scaling_increment <- control$scaling_increment
  pre_diminish <- control$pre_diminish
  
  scale_inc <- scaling_increment / sqrt(max(1, iteration - pre_diminish))
  
  if (any(accept)) {
    scaling[accept] <-
      (sqrt(scaling[accept]) +
         (1 - acceptance_target) * scale_inc) ^ 2
  }
  
  if (any(reject)) {
    scaling[reject] <-
      pmax(sqrt(scaling[reject]) -
             acceptance_target * scale_inc,
           scaling_increment) ^ 2
  }
  scaling
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

update_included <- function(included, i, is_replacement) {
  if (is_replacement) {
    included <- c(included[-1L], i)
  } else {
    included <- c(included, i)
  }
  included
}
