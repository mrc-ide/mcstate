##' Run a pmcmc sampler
##'
##' @title Run a pmcmc sampler
##'
##' @param mcmc_range \code{data.frame} detailing parameters to sample.
##' Must contain columns 'name' (parameter names), 'init' (initial values),
##' 'min' (minimum values), 'max' (maximum values),
##' 'discrete (boolean indicating whether a discrete quantity)' and
##' 'target'
##'
##' @param lprior_funcs functions to calculate log prior for each parameter.
##' A named list for each parameter listed in \code{pars_to_sample}.
##' Each value must be a function which takes named parameter vector as
##'   input, returns a single numeric which is log of the prior probability.
##'
##' @param filter A particle filter to sample
##'
##' @param n_particles Number of particles
##'
##' @param n_mcmc Number of MCMC steps
##'
##' @param proposal_kernel named matrix of proposal covariance for parameters
##'
##' @param run_params List of parameters for \code{particle_filter$run}
##'
##' @param output_proposals Logical indicating whether proposed parameter
##' jumps should be output along with results
##'
##' @param n_chains Number of chains to run
##'
##' @return an mcmc object containing
##' - List of inputs
##' - Matrix of accepted parameter samples, rows = iterations
##'   as well as log prior, (particle filter estimate of) log likelihood
##'   and log posterior
##'
##' @description This is a basic Metropolis-Hastings MCMC sampler.
##' The \code{filter} is run with a set of parameters to evaluate the
##' likeihood. A new set of parameters is proposed, and these likelihoods are
##' compared, jumping with probability equal to their ratio. This is repeated
##' for \code{n_mcmc} proposals, \code{n_chains} independent times.
##'
##' @export
##' @import coda
##' @import furrr
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

  if (!all(par_names %in% names(lprior_funcs))) {
    stop("All sampled parameters must have a defined prior")
  }

  if (!(setequal(rownames(proposal_kernel),
                colnames(proposal_kernel)) &&
         setequal(rownames(proposal_kernel),
                  par_names) &&
      nrow(proposal_kernel) == length(par_names) &&
      ncol(proposal_kernel) == length(par_names))) {
    stop(paste0("proposal_kernel must be a matrix or vector with names ",
                "corresponding to the parameters being sampled"))
  }

  # For the chains to be independent on threads, dust will
  # need to support a long_jump call
  # (for now, could use seed = seed + thread_idx)
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
    names(chains) <- paste0("chain", seq_len(n_chains))

    # calculating rhat
    # convert parallel chains to a coda-friendly format
    chains_coda <- lapply(chains, function(x) {
        traces <- x$results
      coda::as.mcmc(traces[, vars$range$name])
    })

    rhat <- tryCatch(expr = {
      x <- coda::gelman.diag(chains_coda)
      x
    }, error = function(e) {
      print("unable to calculate rhat")
      })


    res <- list(rhat = rhat,
                chains = chains,
                vars = vars)

    class(res) <- "mcstate_pmcmc_list"
  } else {
    res <- chains[[1]]
    class(res) <- "mcstate_pmcmc"
  }

  res
}

# Run a single pMCMC chain (can be done in parallel)
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
    stop("lprior_funcs must return a single numeric representing the log prior")
  }
  if (is.infinite(curr_lprior)) {
    stop("initial parameters are not compatible with supplied prior")
  }

  #
  # Create objects to store outputs
  #

  # initialise output arrays
  res_init <- c(curr_pars,
                "log_prior" = curr_lprior,
                "log_likelihood" = curr_ll,
                "log_posterior" = curr_lpost)
  res <- matrix(data = NA,
                nrow = n_mcmc + 1L,
                ncol = length(res_init),
                dimnames = list(NULL,
                                names(res_init)))
  res[1, ] <- res_init

  if (output_proposals) {
    proposals <- matrix(data = NA,
                        nrow = n_mcmc + 1L,
                        ncol = length(res_init) + 1L,
                        dimnames = list(NULL,
                                        c(names(res_init),
                                          "accept_prob")))
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

  out <- list("results" = as.data.frame(res),
              "acceptance_rate" = 1 - rejection_rate,
              "ess" = ess)

  if (output_proposals) {
    proposals <- as.data.frame(proposals)
    out$proposals <- proposals
  }

 class(out) <- "mcstate_pmcmc"
 out
}

calc_lprior <- function(pars, pars_lprior) {
  lprior <- 0
  for (par in names(pars)) {
      lprior <- lprior + pars_lprior[[par]](pars)
  }
  names(lprior) <- NULL
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

# Check input data.frame
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

  if (range$init < range$min || range$init > range$max) {
    stop("initial parameters are outside of specified range")
  }

  list(range = range,
       index = index)
}

##' @title create a master chain from a pmcmc_list object
##' @param x a pmcmc_list object
##' @param burn_in an integer denoting the number of samples to discard
##' from each chain
##' @export
##'
create_master_chain <- function(x, burn_in) {

  if (class(x) != "mcstate_pmcmc_list") {
    stop("x must be a pmcmc_list object")
  }
  if (!is.numeric(burn_in)) {
    stop("burn_in must be an integer")
  }
  if (burn_in < 0) {
    stop("burn_in must not be negative")
  }
  if (burn_in >= nrow(x$chains[[1]]$results)) {
    stop("burn_in is greater than chain length")
  }

  chains <- lapply(
    X = x$chains,
    FUN = function(z) z$results[-seq_len(burn_in), ]
  )

  master <- do.call(what = rbind, args = chains)
  master
}

#
# Summary and plotting functions
#

##' @export
##' @importFrom stats cor sd
summary.mcstate_pmcmc <- function(object, ...) {

  par_names <- colnames(object$results)
  traces <- object$results

  # calculate correlation matrix
  corr_mat <- round(cor(traces),2)

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

  out <- list("summary" = summ,
              "corr_mat" = corr_mat,
              "sd" = sds)
  out
}

##' @export
summary.mcstate_pmcmc_list <- function(object, ..., burn_in = 1) {

  master_chain <- create_master_chain(x = object,
                                      burn_in = burn_in)

  z <- list(results = master_chain)
  summary.mcstate_pmcmc(z)
}

print_summ <- function(par_name, summ) {
  x <- summ$summary
  paste0(x["mean", par_name],
          "\n(",
          x["2.5%", par_name],
          ", ",
          x["97.5%", par_name], ")")
}

##' @export
##' @importFrom viridis cividis
##' @importFrom graphics hist par plot.new text
plot.mcstate_pmcmc <- function(x, ...) {

  summ <- summary(x)
  par_names <- colnames(x$results)
  n_pars <- length(par_names)

  traces <- x$results
  cols <- viridis::cividis(nrow(traces))
  cols <- cols[order(order(x$results$log_likelihood))]

  par(bty = "n",
      mfcol = c(n_pars, n_pars + 1L),
      mar = c(2.5, 2.5, 2, 1.5),
      mgp = c(1.5, 0.5, 0),
      oma = c(1, 1, 1, 1))

  for (i in seq_len(n_pars)) {
    for(j in seq_len(n_pars)) {
      if (i == j) {
        # plot hists on diagonal
        par_name <- par_names[i]
        breaks <- 10
        hist(traces[[i]],
             main = print_summ(par_name, summ),
             xlab = par_name,
             breaks = breaks,
             cex.main = 1,
             font.main = 1,
             freq = FALSE)
      } else if (i < j) {
        # plot correlations on lower triangle
        plot(x = traces[[i]],
             y = traces[[j]],
             xlab = par_names[i],
             ylab = par_names[j],
             col = cols,
             pch = 20)
      } else if (i > j) {
        # print rho on upper triangle
        plot.new()
        text(x = 0.5,
             y = 0.5,
             labels = paste("r =",
                            summ$corr_mat[i, j]))
      }
    }
  }

  # print traces in final column
  mapply(FUN = plot, traces,
         type = "l",
         ylab = par_names,
         xlab = "Iteration")
}

##' @export
##' @importFrom viridis cividis
##' @importFrom graphics hist par plot.new text lines legend
plot.mcstate_pmcmc_list <- function(x, burn_in = 1, ...) {

  summ <- summary(x, burn_in = burn_in)
  par_names <- colnames(x$chains[[1]]$results)
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

  breaks <- lapply(par_names, function(par_name) {
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
    chain_max <- sapply(h, function(chain) max(chain$density))
    upper_lim <- max(chain_max)
    if (is.na(upper_lim)) {
      upper_lim <- 0
    }
    c(0, upper_lim)
  })

  par(bty = "n",
      mfcol = c(n_pars, n_pars + 1L),
      mar = c(2.5,2.5,1.5,0),
      mgp = c(1.5, 0.5, 0),
      oma = c(1,1,1,1))

  for (i in seq_len(n_pars)) {
    for(j in seq_len(n_pars)) {
      if (i == j) {
        # plot hists on diagonal
        par_name <- par_names[i]
        bs <- breaks[[par_name]]
        plot(x = bs[1],
             y = 1,
             type = "n",
             xlim = c(bs[1], bs[length(bs)]),
             ylim = hist_ylim[[par_name]],
             xlab = par_name,
             ylab = "",
             main = print_summ(par_name, summ),
             cex.main = 1,
             font.main = 1
        )

        mapply(FUN = function(h, col) {
          plot_hists(h = h,
                     col = col,
                     breaks = bs)},
               h = hists[[par_name]],
               col = cols_trace)

      } else if (i < j) {
        # plot correlations on lower triangle
        plot(x = master_chain[[i]],
             y = master_chain[[j]],
             xlab = par_names[i],
             ylab = par_names[j],
             col = cols,
             pch = 20)
      } else if (i > j) {
        # print rho on upper triangle
        plot.new()
        text(x = 0.5,
             y = 0.5, cex = 1.5,
             labels = paste("r =",
                            summ$corr_mat[i, j]))
      }
    }
  }

  # print traces in final column
  n_iter <- nrow(master_chain) / n_chains

  mapply(FUN = function(par_name, leg) {

    trace <- do.call(cbind, traces[[par_name]])
    matplot(x = seq_len(nrow(trace)),
            y = trace,
            type = "l",
            col = cols_trace,
            lty = 1,
            xlab = "Iteration",
            ylab = par_name, )

    if (leg) {
      legend("top",
             ncol = n_chains,
             legend = paste("Chain", seq_len(n_chains)),
             fill = cols_trace,
             bty = "n")
    }
  },
  par_name = par_names,
  leg = c(TRUE, rep(FALSE, length(par_names) - 1)))
}

plot_hists <- function(h, col, breaks) {
  with(h, lines(x =  breaks,
                y = c(density,
                      density[length(density)]),
                type = "s",
                col = col))
}

plot_traces <- function(trace, col) {
  lines(x = seq_along(trace),
        y = trace,
        col = col)
}
