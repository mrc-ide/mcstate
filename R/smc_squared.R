##' Run an SMC^2 algorithm
##'
##' @title Run an SMC^2 algorithm
##'
##' @param data Data to fit to.
##'
##' @param sircovid_model An odin model generator and comparison function
##'
##' @param odin_params Default set of parameters used in odin_model
##'
##' @param beta_params List of default parameters used for time-varying beta
##'
##' @param pars_obs Observation parameters
##'
##' @param pars_seeding Seeding parameters
##'
##' @param earliest_seeding_date Earliest possible date for seeding
##'
##' @param n_particles Number of particles
##'
##' @param n_param_sets Number of parameter sets
##'
##' @param fitted_params list of which parameters are being fitted,
##' with for each a list of log-density function for prior ($dprior)
##' and random generator for prior ($rprior)
##'
##' @param degeneracy_threshold A number in (0,1), threshold used to
##' determine whether or not to do an MCMC move for the parameters
##'
##' @param covariance_scaling Scaling of the weighted covariance for
##' multivariate normal proposal of parameters
##'
##' @export
smc_squared <- function(mcmc_range,
                        lprior_funcs,
                        filter,
                        earliest_seeding_date,
                        n_particles,
                        n_param_sets,
                        fitted_params,
                        degeneracy_threshold = 0.5,
                        covariance_scaling = 0.5) {

  if (n_particles < 2) {
    stop("At least two particles required")
  }
  if (n_param_sets < 2) {
    stop("At least two parameter sets required")
  }

## example:
lprior_funcs <- list(
  beta_start = list(dprior = function(x){dgamma(x,shape=36,scale=3.54e-3,log=TRUE)},
                    rprior = function(){rgamma(1,shape=36,scale=3.54e-3)}
                    ),
  beta_end = list(dprior = function(x){dgamma(x,shape=36,scale=3.54e-3*0.27,log=TRUE)},
                  rprior = function(){rgamma(1,shape=36,scale=3.54e-3*0.27)}
                  )
)

  vars <- mcmc_validate_range(mcmc_range)
  par_names <- as.character(vars$range$name)

  # ALGORITHM
  # Sample n_param_sets of params using rprior
  # Calculate log_prior for these (vectorised version of pMCMC function)
  # Create n_param_sets of filters
  #   [call run2 and set n_steps = 0 (or 1?)]
  # Step each one forward
  #   [call continue with n_steps = 1]
  #   Calculate vector of log likelihoods for each filter
  #   While ESS < covar threshold
  #     Resample these parameter sets
  #     Propose new parameter sets for each filter (vectorised)
  #     Re-run each filter up to the current step with new parameters
  #     Use MH to decide which proposals to accept


  # Intialise fitted parameter sets
  # (these are sampled from their prior distributions)
  params <- initial_params(lprior_funcs, n_param_sets)

  # Calculate log_prior for all parameter sets in params
  log_prior_vec <- calc_log_prior(params, fitted_params)

  # Create filters
  filter_vec <- rep(NA, n_param_sets)
  state <- array(NA_real_, dim = c(length(filter$state(), n_param_sets))
  for (i in n_param_sets) {
    filter_vec[i] <- filter
    filter_vec[i]$run2(n_particles, n_steps = 0,
        save_history = FALSE, vars$index, params)
    state[, i] <- filter_vec[i]$state
  }

  #state <- array(unlist(smc_sq_pre_data_run), 
  #dim = c(length(filter$state()),n_particles,n_param_sets))

  #initialise log likelihood and log weights for each parameter set
  log_likelihood_vec <- rep(0, n_param_sets)
  log_param_weights_vec <- rep(0, n_param_sets)

  n_steps <- filter$n_steps  #TODO extract this from private
  ESS_t <- rep(NA_real_, n_steps)
  acceptance_rate <- rep(NA_real_, n_steps)

  pb <- txtProgressBar(min = 0, max = n_steps, style = 3)
  for (t in seq_len(n_steps) {
    setTxtProgressBar(pb, t-1)
    for (i in n_param_sets) {
      filter_vec[i]$continue(state[, i], pars_compare)
      log_likelihood[i] <- filter_vec[i]$log_likelihood
    }

    #state <- array(unlist(furrr::future_map(model_run, ~ .$state)), dim = dim(state))
    #update log likelihood and log weights for the parameter sets
    #log_ave_weight <- furrr::future_map_dbl(model_run, ~ .$log_ave_weight)
    #log_likelihood <- log_likelihood + log_ave_weight
    #log_param_weights <- log_param_weights + log_ave_weight

    #param_weights are the weights for each parameter set
    #param_weights<-scale_log_weights(log_param_weights)

    param_weights <- scale_log_weights(log_likelihood)
    if (param_weights$average == -Inf) {
      ## Everything is impossible, so stop here
      break
    }

    #Calculate Effective Sample Size
    # (squared-sum of weights/sum of squared-weights)
    ESS <- sum(param_weights$weights)^2/sum(param_weights$weights^2)

    # Resample and MCMC move parameters only if ESS is below a threshold
    if (ESS < degeneracy_threshold * n_param_sets) {

      #calculate weighted covariance of the fitted parameters
      weighted_covar <- cov.wt(params, wt = param_weights$weights)

      #resample from the parameter sets according to their weights
      kappa <- resample(param_weights$weights, "systematic")
      params <- params[kappa, ]
      state <- state[, kappa]
      log_likelihood <- log_likelihood[kappa]
      log_prior <- log_prior[kappa]

      ##MCMC move
      #propose new parameters
      prop_params <- propose_new_params(params,
                      covariance_scaling*weighted_covar$cov)
      #Calculate log prior for the proposed parameter sets
      prop_log_prior <- calc_log_prior(prop_params,fitted_params)

      #set up for below
      prop_log_likelihood <- rep(0,n_param_sets)
      prop_state <- array(NA_real_,dim=dim(state))

      #run the loop below only for possible parameter sets
      kappa <- which(prop_log_prior != -Inf)

      for (k in kappa) {
        filter_vec[k]$run2(n_particles, n_steps = t,
                           save_history = FALSE, vars$index, prop_params)
        log_likelihood[k] <- filter_vec[k]$log_likelihood
      }

      prop_log_likelihood[kappa] <- furrr::future_map_dbl(pf_run, ~ .$log_likelihood)
      prop_state[,,kappa] <- array(unlist(furrr::future_map(pf_run, ~ .$states)), dim = dim(state[,,kappa]))


      #calculate log acceptance probability
      log_accept_prob <- prop_log_likelihood - log_likelihood + prop_log_prior - log_prior
      #which proposed param sets do we accept
      accepted <- which(runif(n_param_sets) < log_accept_prob)
      #update accepted
      params[accepted,] <- prop_params[accepted,]
      log_prior[accepted] <- prop_log_prior[accepted]
      log_likelihood[accepted] <- prop_log_likelihood[accepted]
      state[,,accepted] <- prop_state[,,accepted]

      #reset weights for parameter sets (note that we do not reset when ESS is above threshold)
      log_param_weights[] <- 0
      acceptance_rate[t] <- length(accepted)/n_param_sets
    }
    ESS_t[t]<-ESS
  }

  setTxtProgressBar(pb,t)
  close(pb)

  inputs <- list(data = data,
                 sircovid_model = sircovid_model,
                 odin_params = odin_params,
                 beta_params = beta_params,
                 pars_obs = pars_obs,
                 pars_seeding = pars_seeding,
                 earliest_seeding_date = earliest_seeding_date,
                 n_particles = n_particles,
                 n_param_sets = n_param_sets,
                 fitted_params = fitted_params,
                 degeneracy_threshold = degeneracy_threshold,
                 covariance_scaling = covariance_scaling)

  diagnostics = list(ESS = ESS_t,
                     acceptance_rate = acceptance_rate)

  list(params = params,
       inputs = inputs,
       diagnostics = diagnostics
       )
}


#this function outputs a list of the various parameter objects
#with fitted parameters replaced by values in params
setup_parameters <- function(params, odin_params, beta_params,
                             pars_obs, pars_seeding,
                             generate_beta_func, start_date){

  for (i in seq_len(length(params))){
    if (names(params)[i] %in% names(odin_params)){
      odin_params[[names(params)[i]]] <- params[[names(params)[i]]]
    }
    if (names(params)[i] %in% names(beta_params)){
      beta_params[[names(params)[i]]] <- params[[names(params)[i]]]
    }
    if (names(params)[i] %in% names(pars_obs)){
      pars_obs[[names(params)[i]]] <- params[[names(params)[i]]]
    }
    if (names(params)[i] %in% names(pars_seeding)){
      pars_seeding[[names(params)[i]]] <- params[[names(params)[i]]]
    }
  }

  if ('beta_end' %in% names(beta_params)){
    new_beta <- generate_beta_func(beta_start = beta_params$beta_start,
                                   start_date = start_date,
                                   beta_end = beta_params$beta_end)
  } else if ('beta_reduction' %in% names(beta_params)){
    new_beta <- generate_beta_func(beta_start = beta_params$beta_start,
                                   start_date = start_date,
                                   beta_reduction = beta_params$beta_reduction)
  }

  odin_params$beta_t <- normalise_beta(new_beta$beta_times, odin_params$dt)
  odin_params$beta_y <- new_beta$beta

  list(odin_params = odin_params,
       pars_obs = pars_obs,
       pars_seeding = pars_seeding)
}


## Function to calculate the log_prior given fitted parameters
## pars and (log-)density functions dprior
calc_lprior <- function(pars, lprior_funcs) {
  lprior <- 0
  for (par in names(pars)) {
      lprior <- lprior + lprior_funcs[[par]]$dprior(par)
  }
  names(lprior) <- NULL
  lprior
}


## Initialise fitted parameters by sampling from their prior distributions
initial_params <- function(lprior_funcs, n_param_sets) {

  draw_prior <- function() {
    out <- rep(NA_real_, length(lprior_funcs))
    names(out) <- names(lprior_funcs)
    for (param in lprior_funcs) {
      out[[names(param)[1]]] <- param$rprior()
    }
    unlist(out)
  }

  as.data.frame(t(replicate(n_param_sets, draw_prior(), simplify=TRUE)))
}


#function to propose new parameter sets
##' @importFrom mvtnorm rmvnorm
propose_new_params <- function(params, sigma){

  # propose parameters from a multivariate normal about a mean pars,
  # with covariance matrix sigma
  prop_params <- t(apply(params, 1, FUN = function(pars){rmvnorm(1,mean = as.numeric(pars), sigma = sigma)}))

  #just make sure prop_params is in same form as params
  prop_params <- as.data.frame(prop_params)
  names(prop_params) <- names(params)
  row.names(prop_params) <- NULL

  prop_params

}

#run the model up to the data for the kth set of parameters in params
smc_sq_pre_data_run_model <- function(k, params, data, sircovid_model,
                                odin_params, beta_params, pars_obs,
                                pars_seeding, n_particles,
                                earliest_seeding_date, i_state){

  #Set up the model using the kth set of parameters in params
  pars <- setup_parameters(params[k,], odin_params,
                           beta_params, pars_obs, pars_seeding,
                           sircovid_model$generate_beta_func,
                           earliest_seeding_date)

  odin_model <- sircovid_model$odin_model(user = pars$odin_params)

  seeding_func <- sircovid_model$seeding_model(odin_model,data,pars$pars_seeding)

  X <- pre_data_run(model = odin_model,
                    seeding_func = seeding_func,
                    total_days = NULL,
                    steps_per_day = attr(data, "steps_per_day"),
                    max_seeding_step = data$step_end[[1L]],
                    max_seeding_day = data$day_end[[1L]],
                    n_particles = n_particles,
                    i_state = i_state,
                    save_particles = FALSE)
    X$state
}


#run the model forward and resample for the kth set of parameters in params
smc_sq_run_model_resample <- function(k, t, data, sircovid_model,
                                      state, prev_state, params,
                                      odin_params, beta_params,
                                      pars_obs, pars_seeding,
                                      step, earliest_seeding_date){

  #Set up the model using the kth set of parameters in params
  pars <- setup_parameters(params[k,], odin_params,
                           beta_params, pars_obs, pars_seeding,
                           sircovid_model$generate_beta_func,
                           earliest_seeding_date)

  odin_model <- sircovid_model$odin_model(user = pars$odin_params)

  compare <- sircovid_model$compare_model(odin_model, pars$pars_obs, data)

  state <- particle_run_model(state[,,k], step, odin_model)

  #calculate the log_weights for the kth set, and resample the particles according to those weights
  log_weights <- compare(t, state, prev_state[,,k])
  if (!is.null(log_weights)) {
    weights <- scale_log_weights(log_weights)
    kappa <- resample(weights$weights, "systematic")
    state <- state[,kappa]
    log_ave_weight <- weights$average
  } else{
    log_ave_weight <- 0
  }

  return(list(state = state, log_ave_weight = log_ave_weight))
}


#run a particle filter using the kth set of parameters in prop_params
smc_sq_rerun_pf <- function(k, data, sircovid_model,
                            prop_params, odin_params,
                            beta_params, pars_obs, pars_seeding,
                            n_particles, earliest_seeding_date){

  #Set up the model using the kth set of parameters in params
  pars <- setup_parameters(prop_params[k,], odin_params,
                           beta_params, pars_obs, pars_seeding,
                           sircovid_model$generate_beta_func,
                           earliest_seeding_date)

  odin_model <- sircovid_model$odin_model(user = pars$odin_params)

  compare <- sircovid_model$compare_model(odin_model, pars$pars_obs, data)

  seeding_func <- sircovid_model$seeding_model(odin_model, data, pars$pars_seeding)

  #Run the particle filter for ith proposed parameter set, up until current time point t
  particle_filter(data, odin_model, compare, seeding_func, n_particles, save_end_states = TRUE)
}