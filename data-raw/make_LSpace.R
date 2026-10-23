# Build the attainable (L-CV, L-skewness) regions BurrXII_LSpace, Dagum_LSpace, GGamma_LSpace and
# ExpWeibull_LSpace over the shape boxes of the shared L-moment solver (.lmom_family() in
# R/fitlm_lmom_solver.R). They back LRatio_check() and decide the "outside_domain" fit status of
# fitlm_burr / _dagum / _gengamma / _expweibull, so they must be rebuilt whenever a box changes.
#
# Method (ported from dev/annexD_validation/lspace_v3.R): the box, in the solver's optimisation
# coordinates ((shape1, shape2); (k, d) for the generalised gamma), is covered by a 120 x 120 log
# grid; every grid cell is mapped through the exact L-moment engine to a quadrilateral in
# (L-CV, tau_3) and the quadrilaterals are unioned with sf. Unlike the image of the box edges this
# is exact up to grid resolution even where the map folds (the Jacobian changes sign near the
# degenerate corners). The outer ring of the union is stored as a one-polygon
# sp::SpatialPolygonsDataFrame, the class LRatio_check() reads.
#
# Writes all four objects to R/sysdata.rda (internal) and the exported copies
# data/BurrXII_LSpace.rda and data/Dagum_LSpace.rda. Takes ~10 min (the GG engine dominates).
#
#   Rscript data-raw/make_LSpace.R

suppressMessages(devtools::load_all(".", quiet = TRUE))
N_GRID <- 120

build_lspace <- function(fam) {
  f <- .lmom_family(fam)
  u1 <- exp(seq(log(f$lower[1]), log(f$upper[1]), length.out = N_GRID))
  u2 <- exp(seq(log(f$lower[2]), log(f$upper[2]), length.out = N_GRID))
  L <- array(NA_real_, c(N_GRID, N_GRID, 2))
  for (i in seq_len(N_GRID)) for (j in seq_len(N_GRID)) {
    s <- .lmom_from_coord(fam, c(u1[i], u2[j]))
    v <- as.numeric(f$lmom(1:3, 1, s[1], s[2]))
    L[i, j, ] <- c(v[2] / v[1], v[3])
  }
  quads <- list()
  for (i in seq_len(N_GRID - 1)) for (j in seq_len(N_GRID - 1)) {
    p <- rbind(L[i, j, ], L[i + 1, j, ], L[i + 1, j + 1, ], L[i, j + 1, ], L[i, j, ])
    if (any(!is.finite(p))) next
    q <- sf::st_polygon(list(p)); if (!sf::st_is_valid(q)) q <- sf::st_make_valid(q)
    if (sf::st_area(q) > 0) quads[[length(quads) + 1]] <- q
  }
  u <- sf::st_union(sf::st_sfc(quads))
  u <- sf::st_buffer(sf::st_buffer(u, 1e-4), -1e-4)               # close sliver gaps between quads
  polys <- sf::st_cast(u, "POLYGON")
  outer <- polys[[which.max(sf::st_area(polys))]]
  outer[[1]][, 1:2]
}

# same layout as the historical objects: one Polygons with ID "0", data column FID = 0, no CRS
to_spdf <- function(coords) {
  sp::SpatialPolygonsDataFrame(
    sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords)), ID = "0"))),
    data = data.frame(FID = 0L, row.names = "0"))
}

lspace <- list()
for (fam in c("burr", "dagum", "gengamma", "expweibull")) {
  t0 <- proc.time()[["elapsed"]]
  lspace[[fam]] <- to_spdf(build_lspace(fam))
  cat(sprintf("%-11s %4d vertices, area %.4f, %.0f s\n", fam,
              nrow(lspace[[fam]]@polygons[[1]]@Polygons[[1]]@coords),
              lspace[[fam]]@polygons[[1]]@Polygons[[1]]@area, proc.time()[["elapsed"]] - t0))
}

BurrXII_LSpace    <- lspace$burr
Dagum_LSpace      <- lspace$dagum
GGamma_LSpace     <- lspace$gengamma
ExpWeibull_LSpace <- lspace$expweibull

usethis::use_data(BurrXII_LSpace, Dagum_LSpace, GGamma_LSpace, ExpWeibull_LSpace,
                  internal = TRUE, overwrite = TRUE, compress = "xz")
save(BurrXII_LSpace, file = "data/BurrXII_LSpace.rda", compress = "xz")
save(Dagum_LSpace, file = "data/Dagum_LSpace.rda", compress = "xz")
cat("saved R/sysdata.rda, data/BurrXII_LSpace.rda, data/Dagum_LSpace.rda\n")
