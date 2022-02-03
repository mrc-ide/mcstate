##' Tune pmcmc proposals by running a series of short chains (without
##' saving trajectories etc) and adapting the proposal kernel to match
##' the scaled observed variance-covariance matrix of the parameters.
##'
##' @title Tune pmcmc proposals
##'
##' @param n_tune_blocks The number of tuning blocks.  Using 0 here is
##'   acceptable, resulting in no tuning steps and `initial` being
##'   returned unmodified.  The resulting `vcv` and `proposal_matrix`
##'   elements will be `NA` everywhere, and `pars` will be unmodified.
##'
##' @param n_tune_steps The number of steps per block
##'
##' @param pars A [`pmcmc_parameters`] object containing information
##'   about parameters (ranges, priors, proposal kernel, translation
##'   functions for use with the particle filter). **Important**: the
##'   new proposal matrix will be set into this object (i.e.,
##'   modifying the original)
##'
##' @param control Base [mcstate::pmcmc_control] object.  We will
##'   override a few things in here (especially the number of steps,
##'   trajectory information and burnin/thinning), but we'll reuse
##'   others (parallel control, number of chains, rerun control)
##'
##' @inheritParams pmcmc
##'
##' @param kernel_scaling The scaling to apply to the variance
##'   covariance matrix.  If `NULL` we use the "optimal" scaling
##'   (`2.38^2 / (number of parameters)`)
##'
##' @param weighting_rate The rate at which to reduce weight on
##'   previous blocks of parameters for calculation of the covariance
##'   matrix. Increasing this value towards `Inf` will use only the
##'   last block, while decreasing it to 0 will use all parameters
##'   observed so far.
##'
##' @return A list with elements
##' * `initial`: new initial parameter values, suitable to pass into `pmcmc`
##' * `vcv`: the final estimated variance-covariance matrix
##' * `proposal_kernel`: the final proposal kernel (scaled appropriately)
##' * `history`: a list of parameters and probabilities, as 4d arrays
##'   (parameter, `n_tune_steps`, `n_chains`, `n_tune_blocks`).
##'
##' @export
pmcmc_tune <- function(n_tune_blocks, n_tune_steps,
                       pars, filter, control, initial = NULL,
                       kernel_scaling = NULL, weighting_rate = 1) {
  assert_is(pars, c("pmcmc_parameters", "pmcmc_parameters_nested"))
  assert_is(filter, c("particle_filter", "particle_deterministic"))
  assert_is(control, "pmcmc_control")

  n_pars <- length(pars$names())
  if (is.null(kernel_scaling)) {
    kernel_scaling <- (2.38^2) / n_pars
  }
  assert_scalar_numeric(kernel_scaling)
  assert_scalar_numeric(weighting_rate)

  ## Some manual patching of the control object is required:
  control$n_steps <- n_tune_steps

  ## Disable any in-run state saving
  control$save_state <- FALSE
  control$save_restart <- NULL
  control$save_trajectories <- FALSE

  ## Disable any in-run filtering
  control$n_steps_retain <- n_tune_steps
  control$n_burnin <- 0
  control$n_steps_every <- 1

  history_pars <-
    array(NA_real_, c(n_pars, n_tune_steps, control$n_chains, n_tune_blocks))
  history_probabilities <-
    array(NA_real_, c(3, n_tune_steps, control$n_chains, n_tune_blocks))

  ## In case we make no steps
  vcv <- matrix(NA_real_, n_pars, n_pars,
                dimnames = list(pars$names(), pars$names()))

  for (i in seq_len(n_tune_blocks)) {
    ## TODO: update to use new progress bars, see #199
    message(sprintf("Running tuning stage %d / %d", i, n_tune_blocks))
    results <- pmcmc(pars, filter, control = control, initial = initial)
    history_pars[, , , i] <- t(results$pars)
    history_probabilities[, , , i] <- t(results$probabilities)

    ## Then we look at our pars over time:
    p <- t(array_flatten(history_pars[, , , seq_len(i), drop = FALSE], 2:4))
    wt <- rep(exp(-seq(i - 1, 0) * weighting_rate),
              each = n_tune_steps * control$n_chains)
    vcv <- stats::cov.wt(p, wt)$cov
    pars$update_proposal(vcv * kernel_scaling)

    ## I think there's a question here about how much you want to keep
    ## the chain identity, vs updating things to better points.  We
    ## might make this a tuneable parameter
    initial <- array_drop(
      history_pars[, control$n_steps, , i, drop = FALSE],
      c(2, 4))
  }

  if (!is.null(initial)) {
    if (is.null(dim(initial))) {
      names(initial) <- pars$names()
    } else {
      rownames(initial) <- pars$names()
    }
  }
  rownames(history_pars) <- pars$names()
  rownames(history_probabilities) <-
    c("log_prior", "log_likelihood", "log_posterior")
  history <- list(pars = history_pars,
                  probabilities = history_probabilities)

  list(initial = initial,
       vcv = vcv,
       proposal_kernel = vcv * kernel_scaling,
       history = history)
}
