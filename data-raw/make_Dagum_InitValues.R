# Regenerate data/Dagum_InitValues.rda at UNIT SCALE for the two-step seed in fitlm_dagum.
# One row per unique (shape1, shape2) pair at scale = 1; the scale-free (lcv, tau_3, tau_4) NN
# plus the data-derived scale make the old full scale grid redundant. Shapes from the existing
# table; unit-scale L-moments from the closed-form lmom_dagum engine.
#
# The fitlm_dagum shape box is shape1 in [1e-4, 500], shape2 in [1e-6, 0.95]. The historical
# table stops at shape2 = 0.5, so it is extended with shape1 in log [1e-4, 500] x shape2 in
# (0.5, 0.95] (the grid of dev/annexD_validation/fit_cores_v3.R, seed_table_v3()). Base rows are
# those with shape2 <= 0.5, so re-running the script is idempotent.
#
#   Rscript data-raw/make_Dagum_InitValues.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

pairs <- unique(Dagum_InitValues[Dagum_InitValues$shape2 <= 0.5, c("shape1", "shape2")])
lg <- function(a, b, n) 10^seq(log10(a), log10(b), length.out = n)
ext <- expand.grid(shape1 = lg(1e-4, 500, 40), shape2 = seq(0.5, 0.95, length.out = 20)[-1])
pairs <- unique(rbind(pairs, ext))
cat(sprintf("unique (shape1, shape2) pairs: %d (extension %d)\n", nrow(pairs), nrow(ext)))

lm <- t(mapply(function(s1, s2) {
  v <- anyFit:::lmom_dagum(1:5, scale = 1, shape1 = s1, shape2 = s2)
  c(lambda_1 = v[["lambda_1"]], lambda_2 = v[["lambda_2"]],
    tau_3 = v[["tau_3"]], tau_4 = v[["tau_4"]])
}, pairs$shape1, pairs$shape2))

Dagum_InitValues <- data.frame(
  lambda_1 = lm[, "lambda_1"], lambda_2 = lm[, "lambda_2"],
  tau_3 = lm[, "tau_3"], tau_4 = lm[, "tau_4"],
  lcv = lm[, "lambda_2"] / lm[, "lambda_1"],
  scale = 1, shape1 = pairs$shape1, shape2 = pairs$shape2
)

ok <- with(Dagum_InitValues,
           is.finite(lambda_1) & is.finite(lambda_2) & is.finite(tau_3) & is.finite(tau_4) &
           lcv > 0 & lcv < 1)
Dagum_InitValues <- Dagum_InitValues[ok, ]
cat(sprintf("valid rows kept: %d / %d\n", nrow(Dagum_InitValues), nrow(pairs)))
cat(sprintf("lambda_1 [%.4g, %.4g]; lcv [%.4f, %.4f]\n",
    min(Dagum_InitValues$lambda_1), max(Dagum_InitValues$lambda_1),
    min(Dagum_InitValues$lcv), max(Dagum_InitValues$lcv)))

save(Dagum_InitValues, file = "data/Dagum_InitValues.rda", compress = "xz")
cat("saved data/Dagum_InitValues.rda\n")
