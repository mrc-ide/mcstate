example_sir <- function() {
  model <- odin::odin_(mcstate_file("example/sir/odin_sir.R"), verbose = FALSE)
  sir <- model()
  y0 <- sir$initial()

  compare <- function(state, output, observed, exp_noise = 1e6) {
    incid_modelled <- output[1, ]
    incid_observed <- observed$incid
    lambda <- incid_modelled + rexp(n = length(incid_modelled), rate = exp_noise)
    dpois(x = incid_observed, lambda = lambda, log = TRUE)
  }

  set.seed(1986)
  y <- sir$run(1:100, y0)
  data <- data.frame(step_start = y[, "step"],
                     step_end = y[, "step"] + 1L,
                     incid = y[, "incid"])

  list(model = model, compare = compare, y0 = y0, data = data)
}
