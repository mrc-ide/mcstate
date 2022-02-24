example_sir <- function() {
  set.seed(42)
  model <- dust::dust_example("sir")
  sir <- model$new(pars = list(), step = 0, n_particles = 1)
  y0 <- sir$state()

  compare <- function(state, observed, pars = NULL) {
    if (is.na(observed$incidence)) {
      return(NULL)
    }
    if (is.null(pars$compare$exp_noise)) {
      exp_noise <- 1e6
    } else {
      exp_noise <- pars$compare$exp_noise
    }
    ## This is on the *filtered* state (i.e., returned by run())
    incidence_modelled <- state[1, , drop = TRUE]
    incidence_observed <- observed$incidence
    lambda <- incidence_modelled +
      rexp(n = length(incidence_modelled), rate = exp_noise)
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
  }

  inv_dt <- 4
  day <- seq(1, 100)
  incidence <- rep(NA, length(day))
  history <- array(NA_real_, c(5, 1, length(day) + 1))
  history[, , 1] <- sir$state()

  for (i in day) {
    state_start <- sir$state()
    state_end <- sir$run(i * inv_dt)
    history[, , i + 1] <- state_end
    incidence[i] <- state_end[5, 1]
  }

  data_raw <- data.frame(day = day, incidence = incidence)
  data <- particle_filter_data(data_raw, "day", 4)
  index <- function(info) {
    list(run = 5L, state = 1:3)
  }

  proposal_kernel <- rbind(c(0.00057, 0.00034), c(0.00034, 0.00026))
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

example_volatility <- function(pars = NULL) {
  pars <- pars %||% list(alpha = 0.91, sigma = 1, gamma = 1, tau = 1)

  set.seed(1) # random for init and obs
  volatility <- dust::dust_example("volatility")
  mod <- volatility$new(pars, 0, 1L, seed = 1L)
  mod$update_state(state = matrix(rnorm(1L, 0, 1L), 1))
  steps <- seq(0, 100, by = 1)
  res <- mod$simulate(steps)
  observed <- res[1, 1, -1] + rnorm(length(steps) - 1, 0, 1)
  data <- data.frame(step = steps[-1], observed = observed)
  data <- particle_filter_data(data, "step", 1)

  compare <- function(state, observed, pars) {
    dnorm(observed$observed, pars$compare$gamma * drop(state),
          pars$compare$tau, log = TRUE)
  }

  kalman_filter <- function(pars, data) {
    alpha <- pars$alpha
    sigma <- pars$sigma
    gamma <- pars$gamma
    tau <- pars$tau
    y <- data$observed

    mu <- 0
    s <- 1
    log_likelihood <- 0

    for (t in seq_along(y)) {
      mu <- alpha * mu
      s <- alpha^2 * s + sigma^2
      m <- gamma * mu

      S <- gamma^2 * s + tau^2 # nolint
      K <- gamma * s / S # nolint

      mu <- mu + K * (y[t] - m)
      s <- s - gamma * K * s

      log_likelihood <- log_likelihood + dnorm(y[t], m, sqrt(S), log = TRUE)
    }

    log_likelihood
  }
  list(pars = pars, data = data, compare = compare,
       model = volatility, kalman_filter = kalman_filter)
}

example_sir_shared <- function() {
  set.seed(1)
  model <- dust::dust_example("sir")
  sir <- model$new(pars = list(list(beta = 0.2, gamma = 0.1),
                               list(beta = 0.3, gamma = 0.1)),
                   step = 0, n_particles = 1, pars_multi = TRUE)
  y0 <- sir$state()

  inv_dt <- 4
  day <- seq(1, 100)
  incidence <- matrix(NA, nrow = 2, ncol = length(day))
  history <- array(NA_real_, c(5, 2, length(day) + 1))
  history[, , 1] <- array(y0, c(5, 2, 1))

  for (i in day) {
    state_start <- sir$state()
    state_end <- sir$run(i * inv_dt)
    history[, , i + 1] <- array(state_end, c(5, 2, 1))
    incidence[, i] <- state_end[5, 1, ]
  }

  data_raw <- apply(incidence, 1,
                    function(x) data.frame(day = day, incidence = x))
  data_raw <- do.call(rbind, data_raw)
  data_raw$populations <- factor(rep(letters[1:2], each = nrow(data_raw) / 2))

  data <- particle_filter_data(data_raw, time = "day", rate = 4,
                               population = "populations")

  index <- function(info) {
    list(run = 5L, state = 1:3)
  }

  proposal_fixed <- matrix(0.00026)
  row.names(proposal_fixed) <- colnames(proposal_fixed) <- "gamma"
  proposal_varied <- matrix(0.00057)
  row.names(proposal_varied) <- colnames(proposal_varied) <- "beta"

  pars <- pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("beta", letters[1:2], c(0.2, 0.3),
                                min = 0, max = 1,
                                prior = function(p) log(1e-10)),
         pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                         prior = function(p) log(1e-10))),
    proposal_fixed = proposal_fixed, proposal_varied = proposal_varied)

  compare <- function(state, observed, pars = NULL) {
    if (is.na(observed$incidence)) {
      return(NULL)
    }
    if (is.null(pars$compare$exp_noise)) {
      exp_noise <- 1e6
    } else {
      exp_noise <- pars$compare$exp_noise
    }
    ## This is on the *filtered* state (i.e., returned by run())
    incidence_modelled <- state[1, , drop = TRUE]
    incidence_observed <- observed$incidence
    lambda <- incidence_modelled +
      rexp(n = length(incidence_modelled), rate = exp_noise)
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
  }

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
                           nested = FALSE,
                           state = function() matrix(1, 2, 10),
                           trajectories = function(i) matrix(1, 2, 10)),
                      class = "particle_filter")

  if (is.null(proposal_kernel)) {
    proposal_kernel <- diag(2) * 0.1
    row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("a", "b")
  }

  pars <- pmcmc_parameters$new(
    list(pmcmc_parameter("a", 0.5, min = -1, max = 1,
                         prior = function(p) dunif(p, log = TRUE)),
         pmcmc_parameter("b", 0.5, min = -1, max = 1,
                         prior = function(p) dunif(p, log = TRUE))),
    proposal = proposal_kernel)

  list(target = target, filter = filter, pars = pars)
}


example_uniform_shared <- function(varied = TRUE, fixed = TRUE,
                                   proposal_varied = NULL,
                                   proposal_fixed = NULL) {

  if (!varied || !fixed) {
    n_par <- 2
  } else {
    n_par <- 4
  }

  target <- function(p, ...) {
    rep(1, 3)
  }


  filter <- structure(list(run = target,
                           nested = TRUE,
                           n_particles = 10),
                      class = "particle_filter")

  pars <- list()
  pops <- paste0("p", 1:3)

  if (fixed) {
    if (is.null(proposal_fixed)) {
      proposal_fixed <- diag(2) * 0.1
      row.names(proposal_fixed) <- colnames(proposal_fixed) <- c("a", "b")
    }
    pars <- c(pars,
              list(
                pmcmc_parameter("a", 0.5, min = 0, max = 1,
                                prior = function(p) dunif(p, log = TRUE)),
                pmcmc_parameter("b", 0.5, min = 0, max = 1,
                                prior = function(p) dunif(p, log = TRUE))
              ))
  }

  if (varied) {
    if (is.null(proposal_varied)) {
      proposal_varied <- diag(2) * 0.1
      row.names(proposal_varied) <- colnames(proposal_varied) <- c("c", "d")
    }
    pars <- c(pars,
              list(
                pmcmc_varied_parameter("c", pops, 0.5, min = 0, max = 1,
                                       prior = function(p) dunif(p,
                                                                 log = TRUE)),
                pmcmc_varied_parameter("d", pops, 0.5, min = 0, max = 1,
                                       prior = function(p) dunif(p, log = TRUE))
              ))
  }

  pars <- pmcmc_parameters_nested$new(pars, proposal_varied, proposal_fixed,
                                      pops)

  list(target = target, filter = filter, pars = pars)
}


example_mvnorm <- function() {
  target <- function(p, ...) {
    mvtnorm::dmvnorm(unlist(p), log = TRUE)
  }

  filter <- structure(list(run = target,
                           n_particles = 10,
                           nested = FALSE,
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


example_mvnorm_shared <- function(varied = TRUE, fixed = TRUE,
                                  proposal_varied = NULL,
                                  proposal_fixed = NULL) {
  target <- function(p, ...) {
    vnapply(p, function(x) mvtnorm::dmvnorm(unlist(x), log = TRUE))
  }
  if (!varied || !fixed) {
    n_par <- 2
  } else {
    n_par <- 4
  }
  filter <- structure(list(run = target,
                           nested = TRUE,
                           n_particles = 10),
                      class = "particle_filter")

  pars <- list()
  pops <- paste0("p", 1:3)

  if (fixed) {
    if (is.null(proposal_fixed)) {
      proposal_fixed <- diag(2)
      row.names(proposal_fixed) <- colnames(proposal_fixed) <- c("a", "b")
    }
    pars <- c(pars,
              list(
                pmcmc_parameter("a", 0, min = -100, max = 100),
                pmcmc_parameter("b", 0, min = -100, max = 100)
              ))
  }

  if (varied) {
    if (is.null(proposal_varied)) {
      proposal_varied <- diag(2)
      row.names(proposal_varied) <- colnames(proposal_varied) <- c("c", "d")
    }
    pars <- c(pars,
              list(
                pmcmc_varied_parameter("c", pops, 0, min = -100, max = 100),
                pmcmc_varied_parameter("d", pops, 0, min = -100, max = 100)
              ))
  }

  pars <- pmcmc_parameters_nested$new(pars, proposal_varied, proposal_fixed,
                                      pops)

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
    control <- pmcmc_control(30, save_state = TRUE, save_trajectories = TRUE,
                             save_restart = 40)
    dat$pmcmc <- pmcmc(dat$pars, p, control = control)
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

    control <- pmcmc_control(30, save_state = TRUE, save_trajectories = TRUE,
                             save_restart = 40)

    dat$results <- list(
      pmcmc(dat$pars, p, control = control),
      pmcmc(dat$pars, p, control = control),
      pmcmc(dat$pars, p, control = control))
    test_cache$example_sir_pmcmc2 <- dat
  }
  test_cache$example_sir_pmcmc2
}


example_sir_nested_pmcmc <- function() {
  if (is.null(test_cache$example_sir_nested_pmcmc)) {
    dat <- example_sir_shared()

    n_particles <- 10
    p <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            dat$index, seed = 1L)
    set.seed(1)

    control <- pmcmc_control(30, save_state = TRUE, save_trajectories = TRUE,
                             save_restart = 40)

    dat$results <- list(
      pmcmc(dat$pars, p, control = control),
      pmcmc(dat$pars, p, control = control),
      pmcmc(dat$pars, p, control = control))
    test_cache$example_sir_nested_pmcmc <- dat
  }
  test_cache$example_sir_nested_pmcmc
}


random_array <- function(dim, named = FALSE) {
  if (named) {
    dn <- lapply(seq_along(dim), function(i)
      paste0(LETTERS[[i]], letters[seq_len(dim[i])]))
    names(dn) <- paste0("d", LETTERS[seq_along(dim)])
  } else {
    dn <- NULL
  }
  array(runif(prod(dim)), dim, dimnames = dn)
}


example_variable <- function() {
  ## A small, very silly, model designed to help work with the
  ## multistage filter.  We have a model we can change the dimensions of
  ## without changing the way that the random number draws will work
  ## because only the first entry will be stochastic.
  model <- odin.dust::odin_dust({
    len <- user(integer = TRUE)
    update(x[1]) <- x[1] + rnorm(0, 0.1)
    update(x[2:len]) <- i + step / 10
    initial(x[]) <- 0
    dim(x) <- len
    config(compare) <- "compare_variable.cpp"
  }, verbose = FALSE)

  data <- particle_filter_data(data.frame(time = 1:50, observed = rnorm(50)),
                               "time", 4)
  ## Nonsense model
  compare <- function(state, observed, pars) {
    dnorm(state - observed$observed, log = TRUE)
  }

  index <- function(info) {
    i <- seq(1, info$len, by = 2L)
    names(i) <- letters[i]
    list(run = 1L, state = i)
  }

  transform_state <- function(y, info_old, info_new) {
    n_old <- info_old$len
    n_new <- info_new$len
    if (n_new > n_old) {
      y <- rbind(y, matrix(0, n_new - n_old, ncol(y)))
    } else {
      y <- y[seq_len(n_new), , drop = FALSE]
    }
    y
  }


  list(model = model, data = data, compare = compare, index = index,
       transform_state = transform_state)
}
