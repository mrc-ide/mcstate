# nolint start
# File: RM_RW_odin_test.R
# Purpose: Adaptation of the Ross-McDonald malaria model to aid development
# of parallelisation of deterministic, continuous-time model parameter estimation.
# Includes delay function and time-varying random walk parameter.

################
# Human States #
################
#Rate of change of state variable of human susceptible
deriv(Sh) <- -foi_h * Sh + r * Ih
#Rate of change of the human infectious compartment
deriv(Ih) <- foi_h * Sh - r * Ih
#Force of infection from vectors to humans
foi_h <- a * bh * Iv


#################
# Vector States #
#################
# Rate of change of the susceptible vector compartment
# beta * init_Sv: scale the interpolated emergence to the initial mosquito population size
deriv(Sv) <- beta * init_Sv - foi_v * Sv - mu * Sv
# Latent vector compartments (mutliple compartments to avoid delay function)
deriv(Ev[1]) <- foi_v * Sv - (nrates/tau) * Ev[i] - mu * Ev[i]
deriv(Ev[2:nrates]) <- (nrates/tau) * Ev[i - 1] - (nrates/tau) * Ev[i] - mu * Ev[i]
# Rate of change of the infectious vector population
deriv(Iv) <- (nrates/tau) * Ev[nrates] - mu * Iv
# Force of infection from humans to vectors
foi_v <- a * bv * Ih

# total number of mosquitoes
V <- Sv + sum(Ev[]) + Iv

##################
# Initial States #
##################
initial(Sh) <- init_Sh
initial(Ih) <- init_Ih

initial(Sv) <- init_Sv
initial(Ev[]) <- 0
initial(Iv) <- init_Iv

##############
# User Input #
##############
init_Sh <- 1 - init_Ih
init_Ih <- user()
init_Sv <- user()
init_Iv <- user()

nrates <- user()
dim(Ev) <- nrates

# Ratio mosquitoes:humans
# M <- user(10)

# Biting rate (bites per human per mosquito)
a <- user(1/3)
# Probability of transmission from vectors to humans
bh <- 0.05
# Probability of transmission from humans to vectors
bv <- 0.05
# Daily probability of vector survival
p <- 0.9
mu <- -log(p)
#p <- exp(-mu)
# Rate of recovery
r <- user(1/100)
# Length in mosquito latency period
tau <- user(12)

beta_volatility <- 0.5

# Vector births/deaths - a piece-wise constant function that changes
# every 30 days based on a random-walk function.
initial(beta) <- -log(p)
update(beta) <- beta * exp(rnorm(0, beta_volatility))
#output(beta) <- beta

#Check model equations
N <- Sh + Ih
output(Host_prev) <- Ih / N
output(Vector_prev) <- Iv / V
output(Ev_sum) <- sum(Ev[])
#output(V) <- V
# nolint end