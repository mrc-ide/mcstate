## Definition of the time-step and output as "time"
steps_per_day <- user(1)
dt <- 1 / steps_per_day
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
## compartments, dealing competing hazards for I -> R and I -> D
n_Ix <- rbinom(I, p_IR * dt)
n_IR <- rbinom(n_Ix, 1 - p_death)
n_ID <- n_Ix - n_IR
n_SI <- rbinom(S, p_SI * dt)

## Total population size
N <- S + I + R

## Initial states:
initial(S) <- S_ini
initial(I) <- I_ini
initial(R) <- 0
initial(D) <- 0

## Convert some cumulative values into incidence:
initial(S_inc) <- 0
initial(D_inc) <- 0
update(S_inc) <- if (step %% steps_per_day == 0) n_SI else S_inc + n_SI
update(D_inc) <- if (step %% steps_per_day == 0) n_ID else D_inc + n_ID

## User defined parameters - default in parentheses:
S_ini <- user(1000)
I_ini <- user(10)
beta <- user(0.2)
gamma <- user(0.1)
p_death <- user(0.05)
