##' Run a pmcmc sampler
##'
##' @title Run a pmcmc sampler
##'
##' @param data Data to fit to.  This must be constructed with
##'   \code{particle_filter_data}
##'
##' @param model_params Model parameters, from a call to
##'   \code{generate_parameters()}.
##'
##' @param sircovid_model An odin model generator and comparison function.

##' @param pars_obs list of parameters to use in comparison
##'   with \code{compare_icu}. Must be a list containing, e.g.
##'   list(phi_general = 0.95,
##'   k_general = 2,
##'   phi_ICU = 0.95,
##'   k_ICU = 2,
##'   phi_death = 926 / 1019,
##'   k_death = 2,
##'   exp_noise = 1e6)
##'
##' @param n_mcmc number of mcmc mcmc iterations to perform
##'
##' @param pars_to_sample \code{data.frame} detailing parameters to sample. Must contain columns
##'   'names' (parameter names), 'init' (initial values), 'min' (minimum values), 'max' (maximum values),
##'   'discrete (boolean indicating whether a discrete quantity)'
##'
##' @param pars_lprior functions to calculate log prior for each parameter. A named list for each parameter
##'   listed in \code{pars_to_sample}. Each value must be a function which takes named parameter vector as
##'   input, returns a single numeric which is log of the prior probability.
##'
##' @param proposal_kernel named matrix of proposal covariance for parameters
##'
##' @param n_particles Number of particles
##'
##' @param steps_per_day Number of steps per day
##'
##' @param output_proposals Logical indicating whether proposed parameter jumps should be output along with results
##'
##' @param n_chains Number of chains to run
##'
##' @return an mcmc object containing
##' - List of inputs
##' - Matrix of accepted parameter samples, rows = iterations
##'   as well as log prior, (particle filter estimate of) log likelihood and log posterior
##'
##' @description The user inputs initial parameter values for beta_start and sample_date
##' The log prior likelihood of these parameters is calculated based on the user-defined
##' prior distributions.
##' The log likelihood of the data given the initial parameters is estimated using a particle filter,
##' which has two functions:
##'      - Firstly, to generate a set of 'n_particles' samples of the model state space,
##'        at time points corresponding to the data, one of which is
##'        selected randomly to serve as the proposed state sequence sample at the final
##'        data time point.
##'      - Secondly, to produce an unbiased estimate of the likelihood of the data given the proposed parameters.
##' The log posterior of the initial parameters given the data is then estimated by adding the log prior and
##' log likelihood estimate.
##'
##' The pMCMC sampler then proceeds as follows, for n_mcmc iterations:
##' At each loop iteration the pMCMC sampler perfsorms three steps:
##'   1. Propose new candidate samples for beta_start, beta_end and start_date based on
##'     the current samples, using the proposal distribution
##'     (currently multivariate Gaussian with user-input covariance matrix (proposal_kernel), and reflecting boundaries defined by pars_min, pars_max)
##'   2. Calculate the log prior of the proposed parameters,
##'      Use the particle filter to estimate log likelihood of the data given the proposed parameters, as described above,
##'      as well as proposing a model state space.
##'      Add the log prior and log likelihood estimate to estimate the log posterior of the proposed parameters given the data.
##'   3. Metropolis-Hastings step: The joint canditate sample (consisting of the proposed parameters
##'      and state space) is then accepted with probability min(1, a), where the acceptance ratio is
##'      simply the ratio of the posterior likelihood of the proposed parameters to the posterior likelihood
##'      of the current parameters. Note that by choosing symmetric proposal distributions by including
##'      reflecting boundaries, we avoid the the need to include the proposal likelihood in the MH ratio.
##'
##'   If the proposed parameters and states are accepted then we update the current parameters and states
##'   to match the proposal, otherwise the previous parameters/states are retained for the next iteration.
##'
##' @export
##' @import coda
##' @importFrom stats rnorm
##' @importFrom mvtnorm rmvnorm
##' @importFrom graphics matplot
pmcmc <- function(mcmc_range,
                  lprior_funcs,
                  filter,
                  n_particles,
                  n_mcmc,
                  proposal_kernel,
                  run_params = NULL,
                  output_proposals = FALSE,
                  n_chains = 1) {
  vars <- mcmc_validate_range(mcmc_range)
  par_names <- as.character(vars$range$name)

  if (length(output_proposals) != 1 || !is.logical(output_proposals)) {
    stop("output_proposals must be either TRUE or FALSE")
  }

  if (!all(names(lprior_funcs) %in% par_names)) {
    stop("All sampled parameters must have a defined prior")
  }

  if (!(setequal(rownames(proposal_kernel),
                colnames(proposal_kernel)) &&
         setequal(rownames(proposal_kernel),
                  par_names) &&
      nrow(proposal_kernel) == length(par_names) &&
      ncol(proposal_kernel) == length(par_names))) {
    stop("proposal_kernel must be a matrix or vector with names corresponding
          to the parameters being sampled")
  }

  chains <- furrr::future_pmap(
      .l =  list(n_mcmc = rep(n_mcmc, n_chains)),
      .f = run_mcmc_chain,
      vars,
      lprior_funcs,
      filter,
      n_particles,
      proposal_kernel,
      run_params = run_params,
      output_proposals = output_proposals,
      .progress = TRUE)

  if (n_chains > 1) {
    names(chains) <- paste0('chain', seq_len(n_chains))

    # calculating rhat
    # convert parallel chains to a coda-friendly format
    browser()
    chains_coda <- lapply(chains, function(x) {
        traces <- x$results
      coda::as.mcmc(traces[, vars$range$name])
    })

    rhat <- tryCatch(expr = {
      x <- coda::gelman.diag(chains_coda)
      x
    }, error = function(e) {
      print('unable to calculate rhat')
      })


    res <- list(rhat = rhat,
                chains = chains)

    class(res) <- 'mcstate_pmcmc_list'
  } else {
    res <- chains[[1]]
    class(res) <- 'mcstate_pmcmc'
  }

  res
}

mcmc_validate_range <- function(range) {
  assert_is(range, "data.frame")
  msg <- setdiff(c("name", "init", "min", "max", "discrete", "target"),
                   names(range))
  if (length(msg) > 0L) {
    stop("Missing columns from 'mcmc_range': ",
           paste(squote(msg), collapse = ", "))
  }

  if (anyDuplicated(range$name)) {
    stop("Duplicate 'name' entries not allowed in 'mcmc_range'")
  }

  targets <- c("step_start", "model_data", "pars_compare")
  err <- setdiff(range$target, targets)
  if (length(err) > 0L) {
    stop(sprintf("Invalid target %s: must be one of %s",
                 paste(squote(err), collapse = ", "),
                 paste(squote(targets), collapse = ", ")))
  }

  index <- lapply(targets, function(t) which(range$target == t))
  names(index) <- targets
  if (length(index$step_start) > 1L) {
    stop("At most one target may be 'step_start'")
  }

  if (!is.logical(range$discrete)) {
    stop("'discrete' entries must be TRUE or FALSE")
  }
  if (!is.numeric(range$init)) {
    stop("'init' entries must be numeric")
  }
  if (!is.numeric(range$min)) {
    stop("'min' entries must be numeric")
  }
  if (!is.numeric(range$max)) {
    stop("'max' entries must be numeric")
  }

  if(any(range$init < range$min | range$init > range$mmax)) {
    stop('initial parameters are outside of specified range')
  }

  list(range = range,
       index = index)
}

# Run a single pMCMC chain
run_mcmc_chain <- function(n_mcmc,
                           vars,
                           lprior_funcs,
                           filter,
                           n_particles,
                           proposal_kernel,
                           run_params = NULL,
                           output_proposals = FALSE) {
  #
  # Set initial state
  #
  curr_pars <- vars$range$init
  names(curr_pars) <- vars$range$name

  ## calculate initial prior
  curr_lprior <- calc_lprior(curr_pars, lprior_funcs)

  # run particle filter on initial parameters
  curr_ll <- filter$run2(n_particles,
                         save_history = FALSE,
                         index = vars$index,
                         pars = curr_pars,
                         run_params = run_params)
  curr_lpost <- curr_lprior + curr_ll

  # checks on log_prior and log_likelihood functions
  if (length(curr_lprior) > 1) {
    stop('log_prior must return a single numeric representing the log prior')
  }
  if (is.infinite(curr_lprior)) {
    stop('initial parameters are not compatible with supplied prior')
  }

  #
  # Create objects to store outputs
  #

  # initialise output arrays
  res_init <- c(curr_pars,
                'log_prior' = curr_lprior,
                'log_likelihood' = curr_ll,
                'log_posterior' = curr_lpost)
  res <- matrix(data = NA,
                nrow = n_mcmc + 1L,
                ncol = length(res_init),
                dimnames = list(NULL,
                                names(res_init)))
  res[1, ] <- res_init

  if(output_proposals) {
    proposals <- matrix(data = NA,
                        nrow = n_mcmc + 1L,
                        ncol = length(res_init) + 1L,
                        dimnames = list(NULL,
                                        c(names(res_init),
                                          'accept_prob')))
  }

  #
  # main pmcmc loop
  #
  for (iter in seq_len(n_mcmc) + 1L) {

    # propose new parameters
    prop_pars <- propose_parameters(curr_pars,
                                    proposal_kernel,
                                    vars$range$discrete,
                                    vars$range$min,
                                    vars$range$max)

    ## calculate proposed prior / lhood / posterior
    prop_lprior <- calc_lprior(prop_pars, lprior_funcs)
    prop_ll <- filter$run2(n_particles,
                           save_history = FALSE,
                           index = vars$index,
                           pars = prop_pars,
                           run_params = run_params)
    prop_lpost <- prop_lprior + prop_ll

    # calculate probability of acceptance
    accept_prob <- exp(prop_lpost - curr_lpost)

    # MH step
    if (runif(1) < accept_prob) {
      # update parameters and calculated likelihoods
      curr_pars <- prop_pars
      curr_lprior <- prop_lprior
      curr_ll <- prop_ll
      curr_lpost <- prop_lpost
    }

    # record results
    res[iter, ] <- c(curr_pars,
                     curr_lprior,
                     curr_ll,
                     curr_lpost)

    if (output_proposals) {
      proposals[iter, ] <- c(prop_pars,
                             prop_lprior,
                             prop_ll,
                             prop_lpost,
                             min(accept_prob, 1))
    }

  }

  res <- as.data.frame(res)

  coda_res <- coda::as.mcmc(res)
  rejection_rate <- coda::rejectionRate(coda_res)
  ess <- coda::effectiveSize(coda_res)

  out <- list('results' = as.data.frame(res),
              'acceptance_rate' = 1-rejection_rate,
              "ess" = ess)

  if (output_proposals) {
    proposals <- as.data.frame(proposals)
    out$proposals <- proposals
  }

 class(out) <- 'mcstate_pmcmc'
 out
}

calc_lprior <- function(pars, pars_lprior) {
  lprior <- 0
  for (par in names(pars)) {
      lprior <- lprior + pars_lprior[[par]](pars)
  }
  lprior
}

propose_parameters <- function(pars, proposal_kernel, pars_discrete,
                               pars_min, pars_max) {

  ## proposed jumps are normal with mean pars and sd as input for parameter
  jumps <- pars + drop(rmvnorm(n = 1, sigma = proposal_kernel))

  # discretise if necessary
  jumps[pars_discrete] <- round(jumps[pars_discrete])
  # reflect proposal if it exceeds upper or lower parameter boundary
  jumps <- reflect_proposal(x = jumps,
                          floor = pars_min,
                          cap = pars_max)
  jumps
}

## create function to reflect proposal boundaries at pars_min and pars_max
# this ensures the proposal is symetrical and we can simplify the MH step
reflect_proposal <- function(x, floor, cap) {
  interval <- cap - floor
  abs((x + interval - floor) %% (2 * interval) - interval) + floor
}

##' @title create a master chain from a pmcmc_list object
##' @param x a pmcmc_list object
##' @param burn_in an integer denoting the number of samples to discard from each chain
##' @export
##'
create_master_chain <- function(x, burn_in) {

  if(class(x) != 'pmcmc_list') {
    stop('x must be a pmcmc_list object')
  }
  if(!is.numeric(burn_in)) {
    stop('burn_in must be an integer')
  }
  if(burn_in < 0) {
    stop('burn_in must not be negative')
  }
  if(burn_in >= x$inputs$n_mcmc) {
    stop('burn_in is greater than chain length')
  }

  chains <- lapply(
    X = x$chains,
    FUN = function(z) z$results[-seq_len(burn_in), ]
  )

  do.call(what = rbind, args = chains)
}


##' @export
##' @importFrom stats cor sd
summary.pmcmc <- function(object, ...) {

  par_names <- names(object$inputs$pars$pars_init)

  ## convert start_date to numeric to calculate stats
  data_start_date <- sircovid_date(object$inputs$data$date[1])
  traces <- object$results[,par_names]
  traces$start_date <- sircovid_date(traces$start_date)

  # calculate correlation matrix
  corr_mat <- round(cor(traces),2)


  # add in reduction in beta parameter
  if('beta_end' %in% par_names) {
    traces <- cbind(traces,
                    beta_red = traces$beta_end/traces$beta_start)
  }

  # compile summary
  summ <- rbind(mean = colMeans(traces),
                apply(traces, MARGIN = 2, quantile, c(0.025, 0.975)),
                min = apply(traces, MARGIN = 2, min),
                max =  apply(traces, MARGIN = 2, max)
  )
  summ <- as.data.frame(summ)
  summ <- round(summ, 3)


  sds <- round(apply(traces, 2, sd), 3)
  # convert start_date back into dates
  summ$start_date <- sircovid_date_as_Date(summ$start_date)

  out <- list('summary' = summ,
              'corr_mat' = corr_mat,
              'sd' = sds)
  out

}

##' @export
summary.pmcmc_list <- function(object, ..., burn_in = 101) {

  master_chain <- create_master_chain(x = object,
                                      burn_in = burn_in)

  z <- list(inputs = object$inputs,
            results = master_chain)
  summary.pmcmc(z)

}


##' @export
##' @importFrom viridis cividis
##' @importFrom graphics hist par plot.new text
plot.pmcmc <- function(x, ...) {

  summ <- summary(x)
  par_names <- names(x$inputs$pars$pars_init)

  traces <- x$results[, par_names]
  cols <- viridis::cividis(nrow(traces))
  cols <- cols[order(order(x$results$log_likelihood))]



  par_name <- 'beta_start'
  print_summ <- function(par_name) {
    x <- summ$summary
    paste0(x['mean', par_name],
           '\n(',
           x['2.5%', par_name],
           ', ',
           x['97.5%', par_name], ')')
  }



  n_pars <- length(par_names)


  par( bty = 'n',
       mfcol = c(n_pars, n_pars + 1L),
       mar = c(2.5,2.5,2,1.5),
       mgp = c(1.5, 0.5, 0),
       oma = c(1,1,1,1))


  for (i in seq_len(n_pars)) {
    for(j in seq_len(n_pars)) {

      if (i == j) { # plot hists on diagonal
        par_name <- par_names[i]
        breaks = ifelse(par_name == 'start_date',
                        yes = seq(as.Date('2019-12-01'),
                                  as.Date(x$inputs$data$date[1]), 7),
                        no = 10)
        hist(traces[[i]],
             main = print_summ(par_name),
             xlab = par_name,
             breaks = breaks,
             cex.main = 1,
             font.main = 1,
             freq = FALSE)
      } else if (i < j) {  # plot correlations on lower triangle
        plot(x = traces[[i]],
             y = traces[[j]],
             xlab = par_names[i],
             ylab = par_names[j],
             col = cols,
             pch = 20)
      } else if (i > j) { # print rho on upper triangle
        plot.new()
        text(x = 0.5,
             y=0.5,
             labels = paste('r =',
                            summ$corr_mat[i, j]))
      }
    }
  }

  # print traces in final column
  mapply(FUN = plot, traces,
         type = 'l',
         ylab = par_names,
         xlab = "Iteration")


}

##' @export
##' @importFrom viridis cividis
##' @importFrom graphics hist par plot.new text lines legend
##'
plot.pmcmc_list <- function(x, burn_in = 1, ...) {

  summ <- summary(x, burn_in = burn_in)
  par_names <- names(x$inputs$pars$pars_init)
  n_pars <- length(par_names)

  chains <- x$chains
  n_chains <- length(chains)
  cols_trace <- rev(viridis::viridis(n_chains))


  # compile master chain and order by log posterior for plotting
  master_chain <- create_master_chain(x, burn_in = burn_in)

  master_chain <- master_chain[order(master_chain$log_posterior), ]
  cols <- viridis::cividis(nrow(master_chain))
  cols <- cols[order(master_chain$log_posterior)]




  traces <- lapply(par_names, FUN = function(par_name) {
    lapply(X = chains,
           FUN = function(z) z$results[-seq_len(burn_in), par_name])
  })
  names(traces) <- par_names

  plot_traces <- function(trace, col) {
    lines(x = seq_along(trace),
          y = trace,
          col = col)
  }



  breaks <- lapply(par_names, function(par_name){
    seq(from = min(master_chain[, par_name]),
        to =  max(master_chain[, par_name]),
        length.out = 20)
  })
  names(breaks) <- par_names

  hists <- lapply(par_names, FUN = function(par_name) {
    lapply(X = traces[[par_name]],
           FUN = hist,
           plot = FALSE,
           breaks = breaks[[par_name]])
  })
  names(hists) <- par_names

  hist_ylim <- lapply(hists, function(h) {
    chain_max <- sapply(h, function(chain) max(chain$density) )
    c(0, max(chain_max))
  })

  plot_hists <- function(h, col, breaks) {
    with(h, lines(x =  breaks,
                  y = c(density,
                        density[length(density)]),
                  type = 's',
                  col = col))
  }


  print_summ <- function(par_name) {
    x <- summ$summary
    paste0(x['mean', par_name],
           '\n(',
           x['2.5%', par_name],
           ', ',
           x['97.5%', par_name], ')')
  }


  par( bty = 'n',
       mfcol = c(n_pars, n_pars + 1L),
       mar = c(2.5,2.5,1.5,0),
       mgp = c(1.5, 0.5, 0),
       oma = c(1,1,1,1))


  for (i in seq_len(n_pars)) {
    for(j in seq_len(n_pars)) {

      if (i == j) { # plot hists on diagonal
        par_name <- par_names[i]
        bs <- breaks[[par_name]]
        plot(x = bs[1] ,  # force date axis where needed
             y = 1,
             type = 'n',
             xlim = c(bs[1], bs[length(bs)]),
             ylim = hist_ylim[[par_name]],
             xlab = par_name,
             ylab = '',
             main = print_summ(par_name),
             cex.main = 1,
             font.main = 1
        )

        mapply(FUN = function(h, col) {
          plot_hists(h = h,
                     col = col,
                     breaks = bs)
        },
               h = hists[[par_name]],
               col = cols_trace)


      } else if (i < j) {  # plot correlations on lower triangle
        plot(x = master_chain[[i]],
             y = master_chain[[j]],
             xlab = par_names[i],
             ylab = par_names[j],
             col = cols,
             pch = 20)
      } else if (i > j) { # print rho on upper triangle
        plot.new()
        text(x = 0.5,
             y=0.5, cex = 1.5,
             labels = paste('r =',
                            summ$corr_mat[i, j]))
      }
    }
  }


  # print traces in final column
  n_iter <- nrow(master_chain) / n_chains

  mapply(FUN = function(par_name, leg) {

    trace <- do.call(cbind, traces[[par_name]])
    if(par_name == "start_date") {
      trace <- lubridate::as_date(trace, origin = "1970-01-01")
    }
    matplot(x = seq_len(nrow(trace)),
            y = trace,
            type = "l",
            col = cols_trace,
            lty = 1,
            xlab = 'Iteration',
            ylab = par_name, )

    if(leg) {
      legend('top',
             ncol = n_chains,
             legend = paste('Chain', seq_len(n_chains)),
             fill = cols_trace,
             bty = 'n')
    }
  },
  par_name = par_names,
  leg = c(TRUE, rep(FALSE, length(par_names) - 1)))

}