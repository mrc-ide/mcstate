update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR

initial(S) <- S0
initial(I) <- I0
initial(R) <- 0

p_SI <- 1 - exp(-beta * I / N)
p_IR <- 1 - exp(-gamma)

n_SI <- rbinom(S, p_SI)
n_IR <- rbinom(I, p_IR)

output(incid) <- n_SI

S0 <- user(1000)
I0 <- user(10)
beta <- user(0.2)
gamma <- user(0.1)

N <- S + I + R