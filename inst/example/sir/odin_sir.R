# nolint start
update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR

initial(S) <- S0
initial(I) <- I0
initial(R) <- 0

# Somewhat tedious block here to compute the number infected in this
# block, the total number infected and the total infected over the
# course of a clock-time step (steps_per_day).
initial(infected) <- 0
update(infected) <- n_SI
update(total_infected) <- total_infected + n_SI
initial(total_infected) <- 0
total_infected_prev <- delay(total_infected, steps_per_day)
output(incidence) <- total_infected - total_infected_prev

p_SI <- 1 - exp(-beta * I / N)
p_IR <- 1 - exp(-gamma)

n_SI <- rbinom(S, p_SI * dt)
n_IR <- rbinom(I, p_IR * dt)

output(incid) <- n_SI
output(day) <- step * dt

steps_per_day <- user(4)
S0 <- user(1000)
I0 <- user(10)
beta <- user(0.2)
gamma <- user(0.1)

dt <- 1 / steps_per_day
N <- S + I + R
# nolint end
