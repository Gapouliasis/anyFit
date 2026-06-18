# Closed-form L-moments for the Burr XII and Dagum (Burr III) distributions.
#
# Both reduce the probability-weighted moments beta_r = integral_0^1 x(F) F^r dF to Beta
# functions and then to L-moments via the standard PWM inverse. They return the SAME vector
# layout as lmom::lmrp(order = order, ...): lambda_order for order <= 2, tau_order =
# lambda_order / lambda_2 for order >= 3; and Inf when the mean (hence every L-moment)
# diverges. Derived and validated in dev/ (REPORT.md): exact to 1e-9..1e-13 against two
# independent quantile-/density-integration oracles, and more accurate than lmom::lmrp,
# which fails in the light-tail corner. Internal helpers used by fitlm_burr / fitlm_dagum.

# PWM -> L-moment general inverse: lambda_r = sum_{j=0}^{r-1} (-1)^(r-1-j) C(r-1,j) C(r-1+j,j) beta_j
.pwm_to_lmom <- function(beta) {
  rmax <- length(beta)
  lambda <- numeric(rmax)
  for (r in 1:rmax) {
    j <- 0:(r - 1)
    lambda[r] <- sum((-1)^(r - 1 - j) * choose(r - 1, j) * choose(r - 1 + j, j) * beta[j + 1])
  }
  lambda
}

#' Closed-form L-moments of the Burr XII distribution
#'
#' @keywords internal
#' @noRd
#
# pburr (PW = 1): F(x) = 1 - (1 + a*b*(x/s)^a)^(-1/(a*b)),  a = shape1, b = shape2, s = scale.
# Standard Burr XII F = 1 - (1 + (x/lam)^c)^(-k) with c = a, k = 1/(a*b), lam = s/(a*b)^(1/a);
# c*k = 1/b, so the mean exists iff shape2 < 1.
#   beta_r = lam*k * sum_{j=0}^r (-1)^j C(r,j) Beta(k*(j+1) - 1/c, 1 + 1/c),   c*k > 1
lmom_burr <- function(order = 1:5, scale, shape1, shape2) {
  a <- shape1; b <- shape2; s <- scale
  c_   <- a
  k    <- 1 / (a * b)
  invc <- 1 / c_

  rmax <- max(order)
  out  <- stats::setNames(rep(NA_real_, length(order)),
                          ifelse(order <= 2, paste0("lambda_", order), paste0("tau_", order)))

  # Domain guard: c*k > 1 (shape2 < 1) is required for the mean to exist; k*(j+1) - 1/c > 0
  # for all j >= 0 follows from the j = 0 condition k - 1/c > 0.
  if (k - invc <= 0) {
    out[] <- Inf
    return(out)
  }

  lam <- s / (a * b)^(1 / a)
  beta_r <- function(r) {
    j <- 0:r
    lam * k * sum((-1)^j * choose(r, j) * exp(lbeta(k * (j + 1) - invc, 1 + invc)))
  }
  beta   <- vapply(0:(rmax - 1), beta_r, numeric(1))
  lambda <- .pwm_to_lmom(beta)

  out[] <- ifelse(order <= 2, lambda[order], lambda[order] / lambda[2])
  out
}

#' Closed-form L-moments of the Dagum (Burr III) distribution
#'
#' @keywords internal
#' @noRd
#
# pdagum (PW = 1): F(x) = (1 + (shape2/shape1)*(x/scale)^(-1/shape2))^(-shape1).
# Standard Burr III F = (1 + (x/lam)^(-c))^(-k) with c = 1/shape2, k = shape1,
# lam = scale*(shape2/shape1)^shape2; the mean exists iff shape2 < 1. A single Beta per order:
#   beta_r = lam*k * Beta(k*(r+1) + 1/c, 1 - 1/c),  c > 1  <=>  shape2 < 1   (1/c = shape2)
lmom_dagum <- function(order = 1:5, scale, shape1, shape2) {
  k    <- shape1
  invc <- shape2                       # 1/c = shape2
  lam  <- scale * (shape2 / shape1)^shape2

  rmax <- max(order)
  out  <- stats::setNames(rep(NA_real_, length(order)),
                          ifelse(order <= 2, paste0("lambda_", order), paste0("tau_", order)))

  # Domain guard: 1 - 1/c > 0 (shape2 < 1). The first Beta argument is always > 0.
  if (1 - invc <= 0) {
    out[] <- Inf
    return(out)
  }

  beta   <- vapply(0:(rmax - 1),
                   function(r) lam * k * exp(lbeta(k * (r + 1) + invc, 1 - invc)),
                   numeric(1))
  lambda <- .pwm_to_lmom(beta)

  out[] <- ifelse(order <= 2, lambda[order], lambda[order] / lambda[2])
  out
}
