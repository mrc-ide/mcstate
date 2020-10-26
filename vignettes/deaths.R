## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR - n_ID
update(R) <- R + n_IR
update(D) <- D + n_ID

## Individual probabilities of transition:
p_SI <- 1 - exp(-beta * I / N) # S to I
p_IR <- 1 - exp(-gamma) # I to R

## Draws from binomial distributions for numbers changing between
## compartments:
n_IR <- rbinom(I, p_IR * dt * (1 - p_death))
n_ID <- rbinom(I, p_IR * dt * p_death)
n_SI <- rbinom(S, p_SI * dt)

## Total population size
N <- S + I + R

## Initial states:
initial(S) <- S_ini
initial(I) <- I_ini
initial(R) <- 0
initial(D) <- 0

## User defined parameters - default in parentheses:
S_ini <- user(1000)
I_ini <- user(10)
beta <- user(0.2)
gamma <- user(0.1)
p_death <- user(0.05)
