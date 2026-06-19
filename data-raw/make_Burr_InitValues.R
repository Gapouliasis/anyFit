# Regenerate data/Burr_InitValues.rda at UNIT SCALE for the two-step seed in fitlm_burr.
# The two-step NN matches shapes on the scale-free (lcv, tau_3, tau_4) and then derives the
# scale from the data, so the full scale grid of the old table is redundant: one row per unique
# (shape1, shape2) pair at scale = 1 suffices. Shapes are taken from the existing table; the
# unit-scale L-moments are recomputed with the closed-form lmom_burr engine.
#
#   Rscript data-raw/make_Burr_InitValues.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

pairs <- unique(Burr_InitValues[, c("shape1", "shape2")])
cat(sprintf("unique (shape1, shape2) pairs: %d\n", nrow(pairs)))

lm <- t(mapply(function(s1, s2) {
  v <- anyFit:::lmom_burr(1:5, scale = 1, shape1 = s1, shape2 = s2)
  c(lambda_1 = v[["lambda_1"]], lambda_2 = v[["lambda_2"]],
    tau_3 = v[["tau_3"]], tau_4 = v[["tau_4"]])
}, pairs$shape1, pairs$shape2))

Burr_InitValues <- data.frame(
  lambda_1 = lm[, "lambda_1"], lambda_2 = lm[, "lambda_2"],
  tau_3 = lm[, "tau_3"], tau_4 = lm[, "tau_4"],
  lcv = lm[, "lambda_2"] / lm[, "lambda_1"],
  scale = 1, shape1 = pairs$shape1, shape2 = pairs$shape2
)

ok <- with(Burr_InitValues,
           is.finite(lambda_1) & is.finite(lambda_2) & is.finite(tau_3) & is.finite(tau_4) &
           lcv > 0 & lcv < 1)
Burr_InitValues <- Burr_InitValues[ok, ]
cat(sprintf("valid rows kept: %d / %d\n", nrow(Burr_InitValues), nrow(pairs)))
cat(sprintf("lambda_1 [%.4g, %.4g]; lcv [%.4f, %.4f]\n",
    min(Burr_InitValues$lambda_1), max(Burr_InitValues$lambda_1),
    min(Burr_InitValues$lcv), max(Burr_InitValues$lcv)))

save(Burr_InitValues, file = "data/Burr_InitValues.rda", compress = "xz")
cat("saved data/Burr_InitValues.rda\n")
