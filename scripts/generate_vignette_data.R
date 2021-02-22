set.seed(0)

##
## SIR model
##

# Set up model
gen_sir <- dust::dust_example("sir")
dt <- 0.25
inv_dt <- 1 / dt
sir <- gen_sir$new(data = list(dt = dt, S_ini = 1000, I_ini = 10,
                               beta = 0.2, gamma = 0.1),
                   step = 0,
                   n_particles = 1L,
                   n_threads = 1L,
                   seed = 2L)

## Run model forward 100 days, using 4 steps a day
##
## Calculate the number of new cases by counting how many susceptibles
## have been depleted between updates
n_steps <- 100
incidence <- rep(NA, length.out = n_steps)
true_history <- array(NA_real_, c(length(sir$info()$vars), 1, n_steps + 1))
true_history[, 1, 1] <- sir$state()
for (t in seq_len(n_steps)) {
  state_start <- sir$state()
  state_end <- sir$run(t * inv_dt)
  true_history[, 1, t + 1] <- state_end
  incidence[t] <- state_start[1, 1] - state_end[1, 1]
}

# Add very simple form of noise
incidence <- incidence +
  round(rnorm(mean = 0, n = length(true_history[3, 1, -1]),
              sd = 0.1 * sqrt(true_history[3, 1, -1])))
incidence <- unlist(lapply(incidence, max, 0))
incidence_data <- data.frame(cases = incidence, day = seq_len(n_steps))

# Save outputs
write.csv(incidence_data, "inst/sir_incidence.csv", row.names = FALSE)
saveRDS(true_history, "vignettes/sir_true_history.rds", version = 2)

##
## Model with death data
##
gen_death <- odin.dust::odin_dust("vignettes/deaths.R")

## Set up model
dt <- 0.25
inv_dt <- 1 / dt
sir_death <- gen_death$new(data = list(dt = dt, S_ini = 1000, I_ini = 10,
                                       beta = 0.2, gamma = 0.1),
                           step = 0,
                           n_particles = 1L,
                           n_threads = 1L,
                           seed = 2L)

## Run model forward 100 days, using 4 steps a day
##
## Calculate the number of new cases by counting how many susceptibles
## have been depleted between updates
n_steps <- 100
incidence <- rep(NA, length.out = n_steps)
deaths <- rep(NA, length.out = n_steps)
true_history_deaths <- array(NA_real_,
                             c(length(sir_death$info()$vars), 1, n_steps + 1))
true_history_deaths[, 1, 1] <- sir_death$state()
for (t in seq_len(n_steps)) {
  state_start <- sir_death$state()
  state_end <- sir_death$run(t * inv_dt)
  true_history_deaths[, 1, t + 1] <- state_end
  incidence[t] <- state_start[2, 1] - state_end[2, 1]
  deaths[t] <- state_end[5, 1] - state_start[5, 1]
}

## Add very simple form of noise
missing <- rbinom(n = n_steps, prob = 0.05, size = 1) == 1
incidence[missing] <- NA
death_data <- data.frame(cases = incidence,
                         deaths = deaths,
                         day = seq_len(n_steps))

## Save outputs
write.csv(death_data, "vignettes/death_data.csv", row.names = FALSE)
saveRDS(true_history_deaths, "vignettes/deaths_true_history.rds", version = 2)
