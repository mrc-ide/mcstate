##' Take a grid search produced by \code{\link{grid_search}} and
##' sample \code{n_sample_pars} from the parameter grid uses based
##' on their probability. For each parameter pair chosen, run particle
##' filter with \code{num_particles}.
##'
##' @title Sample Grid Scan
##'
##' @param x Output of \code{\link{grid_search}}.
##'
##' @param ... Other arguments
##'
##' @param filter A \code{particle_filter} object to run (the same
##' (as the one used to produce \code{scan_results})
##'
##' @param n_sample_pars Number of parameter pairs to be sampled. This will
##'   determine how many trajectories are returned. Integer. Default = 10. This
##'   will determine how many trajectories are returned.
##'
##' @param n_particles Number of particles. Positive Integer. Default = 100
##'
##' @param forecast_steps Number of time steps being forecast. Default = 0
##' 
##' @param burn_in Amount of burn in to remove (if an \code{mcstate_pmcmc_list})
##'
##' @return \code{\link{list}}. First element (trajectories) is a 3
##'   dimensional array of trajectories (time, state, tranjectories). Second
##'   element (param_grid) is the parameters chosen when sampling from the
##'   \code{scan_results} grid and the third dimension (inputs) is a list of
##'   model inputs.
##'
##' @export
##' @importFrom utils tail
forecast <- function(x, ...,
                     filter,
                     n_sample_pars = 10,
                     n_particles = 100,
                     forecast_steps = 0,
                     burn_in = 0) {
  assert_integer(n_sample_pars)
  assert_integer(n_particles)
  assert_integer(forecast_steps)
  UseMethod("forecast", x)
}

##' Forecast for grid search
##' @rdname forecast
##' @method forecast mcstate_scan
forecast.mcstate_scan <- function(x, ...,
                                  filter,
                                  n_sample_pars = 10,
                                  n_particles = 100,
                                  forecast_steps = 0) {

  # sample proportional to probability
  sample_idx <- sample(nrow(x$vars$expanded),
                       size = n_sample_pars,
                       replace = TRUE,
                       prob = x$renorm_mat_ll)
  pairs <- x$vars$expanded[sample_idx, ]

  traces <- lapply(split_df_rows(pairs), run_and_forecast,
                   filter, x$vars$index, n_particles,
                   forecast_steps)

  # combine and return
  res <- list("trajectories" = traces,
              "parameters" = pairs)

  class(res) <- "mcstate_forecast"
  return(res)

}

##' Forecast for MCMC
##' @rdname forecast
##' @method forecast mcstate_pmcmc_list
forecast.mcstate_pmcmc_list <- function(x, ...,
                                        filter,
                                        n_sample_pars = 10,
                                        n_particles = 100,
                                        forecast_steps = 0,
                                        burn_in = 1) {

  # discard burn-in and sample
  par_names <- colnames(x$chains[[1]]$results)
  par_names <- par_names[!(par_names %in% c("log_prior", "log_likelihood", "log_posterior"))]
  chains <- create_master_chain(x, burn_in)
  param_grid <- chains[sample.int(n = nrow(chains),
                                  size = n_sample_pars,
                                  replace = FALSE), par_names]

  traces <- lapply(split_df_rows(param_grid), run_and_forecast,
                   filter, x$vars$index, n_particles,
                   forecast_steps)

  # combine and return
  res <- list("trajectories" = traces,
              "parameters" = param_grid)

  class(res) <- "mcstate_forecast"
  return(res)

}

run_and_forecast <- function(model_params, filter, index, n_particles,
                             forecast_steps) {
  filter$run2(n_particles, save_history = TRUE, index, model_params)
  if (forecast_steps > 0) {
    forward_steps <- seq.int(0, forecast_steps)
    trajectories <- filter$predict(forward_steps, append = TRUE)
  } else {
    trajectories <- filter$history
  }
  trajectories
}

##' @export
plot.mcstate_forecast <- function(x, ..., what = "1", data = NULL,
                                  ylab = NULL, title = NULL, col = 'grey80') {
  partition_index = as.integer(what)
  if (partition_index < 1 || partition_index > dim(x$trajectories[[1]])[1]) {
    stop("'what' must be a valid index for a partition")
  }

  if (is.null(title)) {
    title <- paste0("Partition ", partition_index)
  }
  if (is.null(ylab)) {
    ylab <- paste0("Partition ", partition_index)
  }

  particles <- array(0, dim(x$trajectories[[1]])[-1])
  for (i in seq_len(length(x$trajectories))) {
    particles <- particles + x$trajectories[[i]][partition_index, , ]
  }
  particles <- particles / length(x$trajectories)

  plot_particles(particles, ylab = ylab, title = title, col = col)
  if (!is.null(data)) {
    points(data, pch = 19)
  }
}

##' @importFrom graphics plot points matlines
##' @importFrom stats quantile
plot_particles <- function(particles, ylab, title, col = "#44111144") {
  plot(particles[, 1], type = "n", ylab = ylab, xlab = "Step",
       ylim = range(particles, na.rm = TRUE), main = title)
  ## Individual traces
  matlines(t(particles), col=col, lty = 1)
  ## Quantiles
  quantiles <- t(apply(particles, 2, quantile, c(0.025, 0.5, 0.975)))
  matlines(quantiles, col = "black", lty = "dashed")
}