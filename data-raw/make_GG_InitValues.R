# Regenerate data/GG_InitValues.rda at UNIT SCALE for the two-step seed in fitlm_gengamma.
# The old table carried lambda_1 up to ~2.1e93 (extreme shapes), which collapses any absolute
# L-moment NN. The two-step seed matches shapes on the scale-free (lcv, tau_3, tau_4) and derives
# the scale from the data, so one row per unique (shape1, shape2) pair at scale = 1 suffices.
# Shapes from the existing table; unit-scale L-moments from the quadrature lmom_gengamma engine.
# (Supersedes dev/closedform_lmoments/make_GG_InitValues_unit.R, which wrote a dev-only .rds.)
#
#   Rscript data-raw/make_GG_InitValues.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

pairs <- unique(GG_InitValues[, c("shape1", "shape2")])
cat(sprintf("unique (shape1, shape2) pairs: %d\n", nrow(pairs)))

lm <- t(mapply(function(s1, s2) {
  v <- anyFit:::lmom_gengamma(1:5, scale = 1, shape1 = s1, shape2 = s2)
  c(lambda_1 = v[["lambda_1"]], lambda_2 = v[["lambda_2"]],
    tau_3 = v[["tau_3"]], tau_4 = v[["tau_4"]])
}, pairs$shape1, pairs$shape2))

GG_InitValues <- data.frame(
  lambda_1 = lm[, "lambda_1"], lambda_2 = lm[, "lambda_2"],
  tau_3 = lm[, "tau_3"], tau_4 = lm[, "tau_4"],
  lcv = lm[, "lambda_2"] / lm[, "lambda_1"],
  scale = 1, shape1 = pairs$shape1, shape2 = pairs$shape2
)

ok <- with(GG_InitValues,
           is.finite(lambda_1) & is.finite(lambda_2) & is.finite(tau_3) & is.finite(tau_4) &
           lcv > 0 & lcv < 1)
GG_InitValues <- GG_InitValues[ok, ]
cat(sprintf("valid rows kept: %d / %d\n", nrow(GG_InitValues), nrow(pairs)))
cat(sprintf("lambda_1 [%.4g, %.4g]; lcv [%.4f, %.4f]; tau_3 [%.4f, %.4f]\n",
    min(GG_InitValues$lambda_1), max(GG_InitValues$lambda_1),
    min(GG_InitValues$lcv), max(GG_InitValues$lcv),
    min(GG_InitValues$tau_3), max(GG_InitValues$tau_3)))

save(GG_InitValues, file = "data/GG_InitValues.rda", compress = "xz")
cat("saved data/GG_InitValues.rda\n")
