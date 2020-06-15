## This still takes an age to run, especially with 1,000 particles (I
## guess the best-case is going to be 1000 * 63 * 400 steps:
## 25,200,000 steps
##
## The problematic bits seem to be the sampling and shuffling.  We
## might move that into a full C particle filter?
range <- data.frame(name = c("beta", "gamma"),
                    min = c(0.13, 0.05),
                    max = c(0.25, 0.15),
                    n = c(6, 9),
                    target = "pars_model",
                    stringsAsFactors = FALSE)
vars <- grid_search_validate_range(range)

dat <- example_sir()
p <- particle_filter$new(dat$data, dat$model, dat$compare)

i <- 1L
vars$range$target == "pars_model"
vars$expanded[i, ]
vars$index

state <- dat$y0

res <- vnapply(seq_len(nrow(vars$expanded)), function(i)
  p$run2(state, 100, vars$expanded[i, ], vars$index))

res <- array(res, vars$range$n)

image(vars$variables[[1]],
      vars$variables[[2]],
      exp(res - max(res)),
      xlab = names(vars$variables)[1],
      ylab = names(vars$variables)[2])
