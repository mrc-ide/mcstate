

sir_generator <- odin::odin("odin_sir.R")
sir <- sir_generator()

set.seed(1986)
tmp <- sir$run(step = 1:100)


par(mfrow = c(1, 2), bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
plot(tmp[, "incid"], type = "l", ylab = "incidence", xlab = "day")
cols <- 1:3
comps <- c("S", "I", "R")
matplot(tmp[, comps], type = "l", lty = 1,
        ylab = "number of people",
        xlab = "day")
legend("topright", legend = comps, fill = cols, bty = "n", ncol = 3)



data <- tmp[, c("step", "incid")]
colnames(data) <- c("day", "new_cases")
data
write.csv(data, file = "inst/example/sir/sir_data.csv", row.names = FALSE)

compare_output <- function(model, data, exp_noise = 1e6) {
  lambda <- model + rexp(n = length(model), rate = exp_noise)
  dpois(x = data, lambda = lambda, log = TRUE)
}

set.seed(1)
tmp <- sir$run(100)
compare_output(tmp[, "incid"], data$new_cases)
