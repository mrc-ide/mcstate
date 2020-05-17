compare_output <- function(model, data, exp_noise = 1e6) {
  lambda <- model + rexp(n = length(model), rate = exp_noise)
  dpois(x = data, lambda = lambda, log = TRUE)
}

set.seed(1)
tmp <- sir$run(100)
data <- read.csv("inst/example/sir/sir_data.csv")
compare_output(model = tmp[, "incid"], data = data$new_cases)
