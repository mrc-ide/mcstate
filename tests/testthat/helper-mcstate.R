example_sir <- function() {
  model <- dust::dust_example("sir")
  sir <- model$new(data = list(), step = 0, n_particles = 1)
  y0 <- sir$state()

  compare <- function(state, prev_state, observed, pars = NULL) {
    if (is.na(observed$incidence)) {
      return(NULL)
    }
    exp_noise <- pars$compare$exp_noise %||% 1e6
    ## This is on the *filtered* state (i.e., returned by run())
    incidence_modelled <-
      state[1, , drop = TRUE] - prev_state[1, , drop = TRUE]
    incidence_observed <- observed$incidence
    lambda <- incidence_modelled +
      rexp(n = length(incidence_modelled), rate = exp_noise)
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
  }
  inv_dt <- 4
  day <- seq(1, 100)
  incidence <- rep(NA, length(day))
  history <- array(NA_real_, c(4, 1, length(day) + 1))
  history[, , 1] <- sir$state()

  for (i in day) {
    state_start <- sir$state()
    state_end <- sir$run(i * inv_dt)
    history[, , i + 1] <- state_end
    incidence[i] <- state_end[4, 1] - state_start[4, 1]
  }

  data_raw <- data.frame(day = day, incidence = incidence)
  data <- particle_filter_data(data_raw, "day", 4)
  index <- function(info) {
    list(run = 4L, state = 1:3)
  }

  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  pars <- pmcmc_parameters$new(
    list(pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                         prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal = proposal_kernel)

  ## Avoid warnings about scope capture that are not important here.
  environment(index) <- globalenv()
  environment(compare) <- globalenv()

  list(model = model, compare = compare, y0 = y0,
       data_raw = data_raw, data = data, history = history,
       index = index, pars = pars)
}


example_uniform <- function(proposal_kernel = NULL) {
  target <- function(p, ...) {
    1
  }
  filter <- structure(list(run = target,
                           n_particles = 10,
                           state = function() matrix(1, 2, 10),
                           trajectories = function(i) matrix(1, 2, 10)),
                      class = "particle_filter")

  if (is.null(proposal_kernel)) {
    proposal_kernel <- diag(2) * 0.1
    row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("a", "b")
  }

  pars <- pmcmc_parameters$new(
    list(pmcmc_parameter("a", 0.5, min = 0, max = 1),
         pmcmc_parameter("b", 0.5, min = 0, max = 1)),
    proposal = proposal_kernel)

  list(target = target, filter = filter, pars = pars)
}


example_mvnorm <- function() {
  target <- function(p, ...) {
    mvtnorm::dmvnorm(unlist(p), log = TRUE)
  }

  filter <- structure(list(run = target,
                           n_particles = 10,
                           state = function() matrix(1, 2, 10),
                           trajectories = function(i) matrix(1, 2, 10)),
                      class = "particle_filter")

  proposal_kernel <- diag(2)
  pars <- pmcmc_parameters$new(
    list(pmcmc_parameter("a", 0, min = -100, max = 100),
         pmcmc_parameter("b", 0, min = -100, max = 100)),
    proposal = proposal_kernel)

  list(target = target, filter = filter, pars = pars)
}


## Some form of these will likely go back into the package later
acceptance_rate <- function(chain) {
  ## TODO: this is actually pretty awful internally
  1 - coda::rejectionRate(coda::as.mcmc(chain))
}


effective_size <- function(chain) {
  ## TODO: do we ever want the ess of the probabilities?
  coda::effectiveSize(coda::as.mcmc(chain))
}


test_cache <- new.env()


example_sir_pmcmc <- function() {
  if (is.null(test_cache$example_sir_pmcmc)) {
    dat <- example_sir()

    n_particles <- 100
    p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                             index = dat$index)
    set.seed(1)
    dat$pmcmc <- pmcmc(dat$pars, p, 30, TRUE, TRUE)
    test_cache$example_sir_pmcmc <- dat
  }
  test_cache$example_sir_pmcmc
}


example_sir_pmcmc2 <- function() {
  if (is.null(test_cache$example_sir_pmcmc2)) {
    dat <- example_sir()

    n_particles <- 10
    p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                             index = dat$index)
    set.seed(1)

    dat$results <- list(
      pmcmc(dat$pars, p, 30, TRUE, TRUE),
      pmcmc(dat$pars, p, 30, TRUE, TRUE),
      pmcmc(dat$pars, p, 30, TRUE, TRUE))
    test_cache$example_sir_pmcmc2 <- dat
  }
  test_cache$example_sir_pmcmc2
}
