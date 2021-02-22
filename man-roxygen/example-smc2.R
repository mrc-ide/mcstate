# We use an example from dust which implements an epidemiological SIR
# (Susceptible-Infected-Recovered) model; see vignette("sir_models")
# for more background, as this example follows from the pMCMC example
# there

# The key tuning here is the number of particles per filter and number
# of parameter sets to consider simultaneously. Ordinarily these would
# be set (much) higher with an increase in computing time
n_particles <- 42
n_parameter_sets <- 20

# Basic epidemiological (Susceptible-Infected-Recovered) model
sir <- dust::dust_example("sir")

# Pre-computed incidence data
incidence <- read.csv(system.file("sir_incidence.csv", package = "mcstate"))

# Annotate the data so that it is suitable for the particle filter to use
dat <- mcstate::particle_filter_data(incidence, "day", 4)

# Subset the output during run
index <- function(info) {
  list(run = 5L)
}

# The comparison function, used to compare simulated data with observe
# data, given the above subset
compare <- function(state, observed, pars) {
  exp_noise <- 1e6
  incidence_modelled <- state[1L, , drop = TRUE]
  incidence_observed <- observed$cases
  lambda <- incidence_modelled +
    rexp(n = length(incidence_modelled), rate = exp_noise)
  dpois(x = incidence_observed, lambda = lambda, log = TRUE)
}

# Finally, construct the particle filter:
filter <- mcstate::particle_filter$new(dat, sir, n_particles, compare,
                                       index = index)

# To control the smc2 we need to specify the parameters to consider
pars <- mcstate::smc2_parameters$new(
  list(
    mcstate::smc2_parameter("beta",
                            function(n) runif(n, 0, 1),
                            function(x) dunif(x, 0, 1, log = TRUE),
                            min = 0, max = 1),
    mcstate::smc2_parameter("gamma",
                            function(n) runif(n, 0, 1),
                            function(x) dunif(x, 0, 1, log = TRUE),
                            min = 0, max = 1)))
control <- mcstate::smc2_control(n_parameter_sets, progress = TRUE)

# Then we run the particle filter
res <- mcstate::smc2(pars, filter, control)

# This returns quite a lot of information about the fit, and this will
# change in future versions
res

# Most useful is likely the predict method:
predict(res)
