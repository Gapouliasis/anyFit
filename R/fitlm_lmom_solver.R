# Shared method-of-L-moments solver for the four numerically fitted three-parameter families
# (Burr XII, Dagum, Generalised Gamma, Exponentiated Weibull). Ported from the dev prototype
# dev/annexD_validation/fit_cores_v3.R ("v3"); evidence in dev/optimiser_report/report.md and
# dev/gate_v3_realdata/REPORT.md.
#
#   * The L-ratios are scale-free, so only the two shapes are solved for, from the sample
#     (L-CV, tau_3[, tau_4, ...]); scale = lambda_1 / lambda_1(unit scale, shapes) afterwards.
#     There is therefore no scale box.
#   * The shapes are optimised in log space (they span decades; optim's absolute finite-difference
#     step then becomes a relative one), with L-BFGS-B on the plain sum of squared ratio errors
#     (differentiable at the optimum) and factr = 1e6 (tighter values reach the ~1e-14 noise floor
#     of the quadrature engines).
#   * The generalised gamma is optimised in Stacy coordinates (k = shape1 / shape2, d = shape2):
#     its N = 128 Gauss-Legendre engine is inaccurate for k < 0.01, which is a box edge there.
#   * Start: the nearest row of the unit-scale *_InitValues table in min-max scaled
#     (L-CV, tau_3, tau_4). A fit that stalls is restarted once from the second-nearest row.
#   * Every fit gets a status: "matched", "matched_at_bound", "outside_domain" (the sample L-ratios
#     lie outside the attainable region of the family, *_LSpace, and the fit is the least-squares
#     projection onto it) or "not_converged".

# Tolerance on the L-ratio residual for a fit to count as matched.
.LMOM_TOL <- 1e-3

#' Family Definitions for the Shared L-Moment Solver
#'
#' Returns the L-moment engine, seed table, attainable-region polygon and the
#' L-BFGS-B box (in optimisation coordinates: (shape1, shape2), or (k, d) for
#' the generalised gamma) of a numerically fitted family. Built at call time
#' so that the engines and data objects resolve in the package namespace.
#'
#' @param fam One of \code{"burr"}, \code{"dagum"}, \code{"gengamma"},
#'   \code{"expweibull"}.
#' @return A list with elements \code{lmom}, \code{seeds}, \code{lspace},
#'   \code{lower}, \code{upper}.
#' @noRd
.lmom_family <- function(fam) {
  switch(fam,
    burr = list(lmom = lmom_burr, seeds = Burr_InitValues, lspace = BurrXII_LSpace,
                lower = c(0.5, 0.001), upper = c(1e4, 1)),
    dagum = list(lmom = lmom_dagum, seeds = Dagum_InitValues, lspace = Dagum_LSpace,
                 lower = c(1e-4, 1e-6), upper = c(500, 0.95)),
    gengamma = list(lmom = lmom_gengamma, seeds = GG_InitValues, lspace = GGamma_LSpace,
                    lower = c(0.01, 0.05), upper = c(4e4, 200)),          # (k, d)
    expweibull = list(lmom = lmom_expweibull, seeds = ExpWeibull_InitValues,
                      lspace = ExpWeibull_LSpace, lower = c(0.05, 0.02), upper = c(50, 1e4)),
    stop("unknown family: ", fam))
}

# Optimisation coordinates <-> anyFit (shape1, shape2).
.lmom_to_coord   <- function(fam, s) if (fam == "gengamma") c(s[1] / s[2], s[2]) else s
.lmom_from_coord <- function(fam, u) if (fam == "gengamma") c(u[1] * u[2], u[2]) else u

#' Solve the L-Moment Equations of a Three-Parameter Family
#'
#' Fits (scale, shape1, shape2) of a numerically fitted family to sample
#' L-moments: log-shape L-BFGS-B on the scale-free L-ratios selected by
#' \code{order}, analytic scale, one restart on a stall and a fit status.
#'
#' @param fam Family name (see \code{.lmom_family}).
#' @param sample_LM Sample L-moments as returned by \code{Lmoments::Lcoefs}
#'   (\eqn{\lambda_1, \lambda_2, \tau_3, \tau_4, \ldots}), at least up to
#'   order \code{max(4, order)}.
#' @param order Integer vector of L-moment orders to match. Order 2 matches
#'   the L-CV \eqn{\lambda_2 / \lambda_1}, order \eqn{r \ge 3} matches
#'   \eqn{\tau_r}; \eqn{\lambda_1} always sets the scale. At least two ratios
#'   are required.
#' @return A list with \code{par} (scale, shape1, shape2) and \code{diag}
#'   (\code{status}, \code{resid}, \code{convergence}, \code{evals},
#'   \code{restarts}).
#' @noRd
.fit_lmom_shapes <- function(fam, sample_LM, order) {
  f <- .lmom_family(fam)
  ratio_orders <- sort(unique(order[order >= 2]))
  if (length(ratio_orders) < 2)
    stop("`order` must contain at least two L-moment orders >= 2 (order 2 is the L-CV, ",
         "order r >= 3 is tau_r); lambda_1 always sets the scale.")
  rmax <- max(ratio_orders)
  d <- as.numeric(sample_LM)
  lcv <- d[2] / d[1]
  ratios <- function(v) c(v[2] / v[1], v[3:rmax])[ratio_orders - 1]
  target <- ratios(d)
  square <- length(ratio_orders) == 2       # 2 equations, 2 unknowns: an exact match exists inside

  # Two nearest seeds in min-max scaled (lcv, tau_3, tau_4).
  tab <- f$seeds
  rng <- vapply(tab[c("lcv", "tau_3", "tau_4")], function(z) diff(range(z)), numeric(1))
  dist2 <- ((tab$lcv - lcv) / rng[1])^2 + ((tab$tau_3 - d[3]) / rng[2])^2 +
           ((tab$tau_4 - d[4]) / rng[3])^2
  nn <- order(dist2)[1:2]

  obj <- function(lu) {
    s <- .lmom_from_coord(fam, exp(lu))
    v <- as.numeric(f$lmom(1:rmax, 1, s[1], s[2]))
    if (any(!is.finite(v))) return(1e6)
    e <- sum((ratios(v) - target)^2)
    if (is.finite(e)) e else 1e6
  }
  run <- function(i) {
    u0 <- pmin(pmax(.lmom_to_coord(fam, as.numeric(tab[i, c("shape1", "shape2")])), f$lower),
               f$upper)
    o <- stats::optim(log(u0), obj, method = "L-BFGS-B", lower = log(f$lower),
                      upper = log(f$upper), control = list(factr = 1e6, maxit = 500))
    u <- exp(o$par)
    o$shapes <- .lmom_from_coord(fam, u)
    o$at_bound <- any(abs(u - f$lower) <= 1e-6 * f$lower) || any(abs(u - f$upper) <= 1e-6 * f$upper)
    v <- as.numeric(f$lmom(1:rmax, 1, o$shapes[1], o$shapes[2]))
    o$resid <- max(abs(ratios(v) - target))
    o
  }

  # Square system: the residual is the convergence test and drives the restart. Least squares on
  # 3+ ratios: the residual is a sampling-noise compromise, so only optim's code 52
  # (ABNORMAL_TERMINATION_IN_LNSRCH) restarts and the status follows the optimiser code.
  o <- run(nn[1]); restarts <- 0L; evals <- unname(o$counts[1])
  stalled <- if (square) !isTRUE(o$resid < .LMOM_TOL) && (o$convergence == 52 || !o$at_bound) else
    o$convergence == 52
  if (stalled) {
    o2 <- run(nn[2]); restarts <- 1L; evals <- evals + unname(o2$counts[1])
    better <- if (square) isTRUE(o2$resid < o$resid) else o2$convergence == 0 && o$convergence != 0
    if (better) o <- o2
  }

  s <- o$shapes
  scale <- d[1] / as.numeric(f$lmom(1, 1, s[1], s[2]))[1]
  poly <- f$lspace@polygons[[1]]@Polygons[[1]]@coords
  inside <- isTRUE(sp::point.in.polygon(lcv, d[3], poly[, 1], poly[, 2]) > 0)
  status <- if (square) {
    if (!is.finite(o$resid)) "not_converged" else
    if (o$resid < .LMOM_TOL) { if (o$at_bound) "matched_at_bound" else "matched" } else
    if (!inside) "outside_domain" else "not_converged"
  } else {
    if (o$convergence != 0) "not_converged" else if (!inside) "outside_domain" else
    if (o$at_bound) "matched_at_bound" else "matched"
  }
  list(par = c(scale, s),
       diag = list(status = status, resid = o$resid, convergence = o$convergence,
                   evals = evals, restarts = restarts))
}
