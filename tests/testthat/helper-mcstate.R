models <- list(
  sir = odin::odin_(mcstate_file("example/sir/odin_sir.R"), verbose = FALSE))


example_sir <- function() {
  model <- models$sir
  sir <- model()
  y0 <- sir$initial(0)

  compare <- function(state, output, observed, exp_noise = 1e6) {
    incidence_modelled <- output[1, ]
    incidence_observed <- observed$incidence
    lambda <- incidence_modelled +
      rexp(n = length(incidence_modelled), rate = exp_noise)
    dpois(x = incidence_observed, lambda = lambda, log = TRUE)
  }

  set.seed(1986)
  y <- sir$run(seq(0, 400, by = 4), y0)

  data_raw <- as.data.frame(y)[c("day", "incidence")]
  data <- particle_filter_data(data_raw[-1, ], "day", 4)

  list(model = model, compare = compare, y0 = y0,
       y = y, data_raw = data_raw, data = data)
}
