##' Take a grid search produced by \code{\link{grid_search}} and
##' sample \code{n_sample_pairs} from the parameter grid uses based
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
##' @param n_sample_pairs Number of parameter pairs to be sampled. This will
##'   determine how many trajectories are returned. Integer. Default = 10. This
##'   will determine how many trajectories are returned.
##'
##' @param n_particles Number of particles. Positive Integer. Default = 100
##'
##' @param forecast_steps Number of time steps being forecast. Default = 0
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
                     n_sample_pairs = 10,
                     n_particles = 100,
                     forecast_steps = 0) {
  UseMethod("forecast", x)
}

##' Forecast for grid search
##' @rdname forecast
##' @method forecast mcstate_scan
forecast.mcstate_scan <- function(x, ...,
                                  filter,
                                  n_sample_pairs = 10,
                                  n_particles = 100,
                                  forecast_steps = 0) {

  # checks on args
  assert_integer(n_sample_pairs)
  assert_integer(n_particles)

  # sample proportional to probability
  sample_idx <- sample(nrow(x$vars$expanded),
                       size = n_sample_pairs,
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

##' Take a parameter search produced by \code{\link{pmcmc}} and
##' sample \code{n_sample} from the parameter space
##' on their probability. For each parameter pair chosen, run particle
##' filter with \code{num_particles} and sample 1 trajectory
##'
##' @title Sample pmcmc
##'
##' @param mcmc_results Output of \code{\link{pmcmc}}.
##'
##' @param burn_in Number of burn-in samples to discard
##'
##' @param n_sample Number of parameter pairs to be sampled. This will
##'   determine how many trajectories are returned. Integer. Default = 10. This
##'   will determine how many trajectories are returned.
##'
##' @param n_particles Number of particles. Positive Integer. Default = 100
##'
##' @param forecast_days Number of days being forecast. Default = 0
##'
##' @return \code{\link{list}}. First element (trajectories) is a 3
##'   dimensional array of trajectories (time, state, tranjectories). Second
##'   element (param_grid) is the parameters chosen when  sampling from the
##'   \code{scan_results} grid and the third dimension (inputs) is a list of
##'   model inputs.
##'
##' @export
##' @importFrom furrr future_map
##' @importFrom purrr transpose
forecast.mcstate_pmcmc <- function(mcmc_results,
                         burn_in = 1,
                         n_sample = 10,
                         n_particles = 100,
                         forecast_days = 0) {

  # checks on args
  assert_pos_int(n_sample)
  assert_pos_int(n_particles)
  assert_pos_int(forecast_days)

  # discard burn-in
  pars_to_sample <- names(mcmc_results$inputs$pars$pars_init)
  chains <- create_master_chain(mcmc_results, burn_in)
  param_grid <- chains[sample.int(n = nrow(chains),
                                  size  = n_sample,
                                  replace = FALSE), pars_to_sample]

  if ("beta_changepoints" %in% names(mcmc_results$inputs)){
    beta_changepoints <- mcmc_results$inputs$beta_changepoints
  } else {
    beta_changepoints <- NULL
  }

  forecasts <- function(sampled_pars) {
    pars <- as.list(sampled_pars)
    pars$start_date <- sircovid_date(pars$start_date)
    trace <- calc_loglikelihood(pars,
                                mcmc_results$inputs$data,
                                mcmc_results$inputs$sircovid_model,
                                mcmc_results$inputs$model_params,
                                mcmc_results$inputs$steps_per_day,
                                mcmc_results$inputs$pars_obs,
                                n_particles,
                                forecast_days,
                                return = "full",
                                beta_changepoints)
    trace
  }

  # Run forecasts in parallel
  traces <- furrr::future_map(.x = purrr::transpose(param_grid), .f = forecasts)
  trajectories <- traces_to_trajectories(traces)

  if (is.null(beta_changepoints)){
    # combine and return
    res <- list("trajectories" = trajectories,
                "param_grid" = param_grid,
                inputs = list(
                  model_params = mcmc_results$inputs$model_params,
                  pars_obs = mcmc_results$inputs$pars_obs,
                  data = mcmc_results$inputs$data,
                  model = mcmc_results$inputs$sircovid_model,
                  forecast_days = forecast_days))
  } else {
    res <- list("trajectories" = trajectories,
                "param_grid" = param_grid,
                inputs = list(
                  model_params = mcmc_results$inputs$model_params,
                  pars_obs = mcmc_results$inputs$pars_obs,
                  data = mcmc_results$inputs$data,
                  model = mcmc_results$inputs$sircovid_model,
                  forecast_days = forecast_days,
                  beta_changepoints = beta_changepoints))
  }

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