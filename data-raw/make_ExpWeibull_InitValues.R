# Build data/ExpWeibull_InitValues.rda at UNIT SCALE for the two-step seed in fitlm_expweibull.
# Unlike Burr/Dagum/GG there is no pre-existing ExpWeibull table to take shapes from, so the
# (shape1, shape2) grid is built FRESH over the fitlm_expweibull optim box (shape1 in [0.05,50],
# shape2 in [0.02,50]); the two-step NN matches shapes on the scale-free (lcv, tau_3, tau_4) and
# derives the scale from the data, so one row per (shape1, shape2) at scale = 1 suffices. The
# unit-scale L-moments are computed with the closed-form tanh-sinh engine anyFit:::lmom_expweibull.
#
#   Rscript data-raw/make_ExpWeibull_InitValues.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

# Log-spaced shape grid spanning the optim box; validation-grade fixed rule (h, X).
shape1_grid <- exp(seq(log(0.05), log(50), length.out = 60))   # a = shape1
shape2_grid <- exp(seq(log(0.02), log(50), length.out = 60))   # b = shape2
pairs <- expand.grid(shape1 = shape1_grid, shape2 = shape2_grid)
cat(sprintf("(shape1, shape2) grid pairs: %d\n", nrow(pairs)))

lm <- t(mapply(function(s1, s2) {
  v <- anyFit:::lmom_expweibull(1:5, scale = 1, shape1 = s1, shape2 = s2, h = 1/64, X = 5)
  c(lambda_1 = v[["lambda_1"]], lambda_2 = v[["lambda_2"]],
    tau_3 = v[["tau_3"]], tau_4 = v[["tau_4"]])
}, pairs$shape1, pairs$shape2))

ExpWeibull_InitValues <- data.frame(
  lambda_1 = lm[, "lambda_1"], lambda_2 = lm[, "lambda_2"],
  tau_3 = lm[, "tau_3"], tau_4 = lm[, "tau_4"],
  lcv = lm[, "lambda_2"] / lm[, "lambda_1"],
  scale = 1, shape1 = pairs$shape1, shape2 = pairs$shape2
)

ok <- with(ExpWeibull_InitValues,
           is.finite(lambda_1) & is.finite(lambda_2) & is.finite(tau_3) & is.finite(tau_4) &
           lcv > 0 & lcv < 1)
ExpWeibull_InitValues <- ExpWeibull_InitValues[ok, ]
cat(sprintf("valid rows kept: %d / %d\n", nrow(ExpWeibull_InitValues), nrow(pairs)))
cat(sprintf("lcv [%.4f, %.4f]; tau_3 [%.4f, %.4f]; tau_4 [%.4f, %.4f]\n",
    min(ExpWeibull_InitValues$lcv),   max(ExpWeibull_InitValues$lcv),
    min(ExpWeibull_InitValues$tau_3), max(ExpWeibull_InitValues$tau_3),
    min(ExpWeibull_InitValues$tau_4), max(ExpWeibull_InitValues$tau_4)))

save(ExpWeibull_InitValues, file = "data/ExpWeibull_InitValues.rda", compress = "xz")
cat("saved data/ExpWeibull_InitValues.rda\n")
