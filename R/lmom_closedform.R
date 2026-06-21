# Closed-form L-moments for the Burr XII and Dagum (Burr III) distributions.
#
# Both reduce the probability-weighted moments beta_r = integral_0^1 x(F) F^r dF to Beta
# functions and then to L-moments via the standard PWM inverse. They return the SAME vector
# layout as lmom::lmrp(order = order, ...): lambda_order for order <= 2, tau_order =
# lambda_order / lambda_2 for order >= 3; and Inf when the mean (hence every L-moment)
# diverges. Derived and validated in dev/ (REPORT.md): exact to 1e-9..1e-13 against two
# independent quantile-/density-integration oracles, and more accurate than lmom::lmrp,
# which fails in the light-tail corner. Internal helpers used by fitlm_burr / fitlm_dagum.

#' Convert Probability-Weighted Moments to L-Moments
#'
#' Applies the general PWM-to-L-moment inverse:
#' \eqn{\lambda_r = \sum_{j=0}^{r-1} (-1)^{r-1-j} \binom{r-1}{j} \binom{r-1+j}{j} \beta_j}.
#'
#' @param beta Numeric vector of probability-weighted moments
#'   \eqn{\beta_0, \beta_1, \ldots}.
#'
#' @return Numeric vector of L-moments \eqn{\lambda_1, \lambda_2, \ldots}.
#' @noRd
.pwm_to_lmom <- function(beta) {
  rmax <- length(beta)
  lambda <- numeric(rmax)
  for (r in 1:rmax) {
    j <- 0:(r - 1)
    lambda[r] <- sum((-1)^(r - 1 - j) * choose(r - 1, j) * choose(r - 1 + j, j) * beta[j + 1])
  }
  lambda
}

#' Closed-Form L-Moments of the Burr Type XII Distribution
#'
#' Computes L-moments of the Burr XII distribution via closed-form Beta-function
#' expressions for the probability-weighted moments. The mean and all higher
#' L-moments exist only when \code{shape2 < 1}; the function returns \code{Inf}
#' outside this domain. Exact to machine precision and more accurate than
#' \code{lmom::lmrp} in the light-tail corner.
#'
#' @param order Integer vector of L-moment orders (e.g. \code{1:5}).
#' @param scale Scale parameter (\eqn{s > 0}).
#' @param shape1 First shape parameter (\eqn{\zeta > 0}).
#' @param shape2 Second shape parameter (\eqn{\theta > 0}); the mean exists
#'   only when \code{shape2 < 1}.
#'
#' @return Named numeric vector of L-moments (\code{lambda_1}, \code{lambda_2},
#'   \code{tau_3}, \code{tau_4}, ...) or \code{Inf} when the mean diverges.
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

#' Closed-Form L-Moments of the Dagum (Burr III) Distribution
#'
#' Computes L-moments of the Dagum distribution via closed-form Beta-function
#' expressions for the probability-weighted moments. The mean and all higher
#' L-moments exist only when \code{shape2 < 1}; the function returns \code{Inf}
#' outside this domain.
#'
#' @param order Integer vector of L-moment orders (e.g. \code{1:5}).
#' @param scale Scale parameter (\eqn{s > 0}).
#' @param shape1 First shape parameter (\eqn{\alpha > 0}), controls upper-tail
#'   heaviness.
#' @param shape2 Second shape parameter (\eqn{\theta > 0}); tail index — the mean
#'   exists only when \code{shape2 < 1}.
#'
#' @return Named numeric vector of L-moments (\code{lambda_1}, \code{lambda_2},
#'   \code{tau_3}, \code{tau_4}, ...) or \code{Inf} when the mean diverges.
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

# Gauss-Legendre nodes/weights on [-1,1] via Golub-Welsch (base R; no statmod dep).
# Eigen-decomposition of the symmetric tridiagonal Jacobi matrix; nodes = eigenvalues,
# weights = 2*(first eigenvector component)^2. Cached by N.

#' Cache Environment for Gauss-Legendre Rules
#'
#' Environment storing pre-computed Gauss-Legendre nodes and weights keyed by
#' the number of nodes N. Populated lazily by \code{.gl_rule}.
#'
#' @noRd
.gl_cache <- new.env(parent = emptyenv())

#' Gauss-Legendre Nodes and Weights (Golub-Welsch)
#'
#' Computes N-point Gauss-Legendre quadrature nodes and weights on \eqn{[-1, 1]}
#' via the Golub-Welsch algorithm (eigen-decomposition of the symmetric
#' tridiagonal Jacobi matrix). Results are cached in \code{.gl_cache} by N.
#'
#' @param N Number of quadrature nodes.
#'
#' @return A list with elements \code{x} (nodes on \eqn{[-1, 1]}) and \code{w}
#'   (weights).
#' @noRd
.gl_rule <- function(N) {
  key <- as.character(N)
  if (is.null(.gl_cache[[key]])) {
    i <- seq_len(N - 1)
    bk <- i / sqrt(4 * i^2 - 1)
    J <- matrix(0, N, N)
    J[cbind(i, i + 1)] <- bk
    J[cbind(i + 1, i)] <- bk
    e <- eigen(J, symmetric = TRUE)
    o <- order(e$values)
    .gl_cache[[key]] <- list(x = e$values[o], w = 2 * (e$vectors[1, o])^2)
  }
  .gl_cache[[key]]
}

#' Quadrature L-Moments of the Stacy Generalised Gamma Distribution
#'
#' Computes L-moments of the generalised gamma distribution via Gauss-Legendre
#' quadrature of the probability-weighted moments. A change of variable maps
#' the PWM integrals to the smooth, bounded unit interval, eliminating endpoint
#' singularities and producing a continuous, optimiser-friendly surface.
#' A fixed 128-node Gauss-Legendre rule yields accuracy far below sampling
#' variability.
#'
#' @param order Integer vector of L-moment orders (e.g. \code{1:5}).
#' @param scale Scale parameter (\eqn{s > 0}).
#' @param shape1 First shape parameter (\eqn{\alpha > 0}); the Stacy shape is
#'   \eqn{k = \alpha/\beta}.
#' @param shape2 Second shape parameter (\eqn{\beta > 0}); \eqn{\beta = 1}
#'   recovers the ordinary gamma distribution.
#' @param N Number of Gauss-Legendre nodes. Default 128.
#'
#' @return Named numeric vector of L-moments (\code{lambda_1}, \code{lambda_2},
#'   \code{tau_3}, \code{tau_4}, ...).
#' @keywords internal
#' @noRd
#
# dgengamma/pgengamma/qgengamma use VGAM::*gengamma.stacy(scale, k = shape1/shape2, d = shape2),
# so the sigma=1 quantile is Q1(u) = qgamma(u, shape = k, rate = 1)^(1/d). The quantile is the
# inverse incomplete gamma -- no elementary closed form -- so the L-moments are computed by
# RIGOROUS Gauss-Legendre quadrature of the PWMs beta_j = integral_0^1 Q1(u) u^j du. Evaluating
# in the gamma-latent variable Y = X^d ~ Gamma(k,1) and substituting v = pgamma(y, k+1/d) maps
# beta_j to a BOUNDED, SMOOTH, singularity-free unit-interval integral:
#   beta_j = mean0 * integral_0^1 [ pgamma(qgamma(v, k+1/d), k) ]^j dv ,  mean0 = Gamma(k+1/d)/Gamma(k)
# A fixed N=128 GL rule gives ~1e-6 on lambda1/lambda2 (far below sampling noise) and, crucially,
# makes the L-moments a SMOOTH function of (scale,shape1,shape2) for the finite-difference optim.
# stats:: qualifier is mandatory: the anyFit namespace masks qgamma/pgamma/dgamma.
lmom_gengamma <- function(order = 1:5, scale, shape1, shape2, N = 128) {
  k <- shape1 / shape2
  d <- shape2
  rmax <- max(order)
  out <- stats::setNames(rep(NA_real_, length(order)),
                         ifelse(order <= 2, paste0("lambda_", order), paste0("tau_", order)))

  if (!is.finite(k) || !is.finite(d) || k <= 0 || d <= 0) return(out)

  ph    <- 1 / d
  mean0 <- exp(lgamma(k + ph) - lgamma(k))          # = beta_0 (sigma = 1), exact
  if (!is.finite(mean0)) return(out)

  rseq <- 0:(rmax - 1)
  gl <- .gl_rule(N)
  v  <- 0.5 + 0.5 * gl$x; w <- 0.5 * gl$w
  G  <- stats::pgamma(stats::qgamma(v, shape = k + ph), shape = k)   # bounded in [0,1]
  beta <- mean0 * vapply(rseq, function(j) sum(w * G^j), numeric(1))
  if (!all(is.finite(beta))) return(out)

  lambda <- .pwm_to_lmom(beta) * scale

  out[] <- ifelse(order <= 2, lambda[order], lambda[order] / lambda[2])
  out
}

# ---- Exponentiated (Generalized) Weibull --------------------------------------------------
# pexpweibull: F(x) = [1 - exp(-(x/scale)^shape1)]^shape2,  a = shape1, b = shape2, s = scale.
# Quantile x(F) = s*(-log(1 - F^(1/b)))^(1/a). The PWMs beta_r = integral_0^1 x(F) F^r dF reduce
# to the exact generalised-binomial series beta_r = s*b*Gamma(1+1/a) sum_k (-1)^k C(b(r+1)-1,k) /
# (k+1)^(1+1/a), which is catastrophically ill-conditioned for large b(r+1)-1. Analytically
# resumming it gives the WELL-CONDITIONED integral
#   beta_r = s*b * integral_0^1 y^(b(r+1)-1) * (-log(1-y))^(1/a) dy,
# evaluated with a fixed tanh-sinh (double-exponential) rule whose nodes/weights vanish doubly-
# exponentially at both endpoints -> ~1e-12 for ALL a,b,s>0 with pure O(N) arithmetic, no per-call
# adaptive integration. All L-moments exist for any a,b,s>0, so there is no domain guard. b=1 is the
# ordinary Weibull; a=b=1 the exponential (tau3=1/3, tau4=1/6, tau5=1/10). Derived and validated in
# dev/closedform_lmoments (expweibull_implementation_vs_methodology.md): machine precision vs two
# independent quantile-/density-integration oracles. Used by fitlm_expweibull.

# tanh-sinh nodes/weights for integral_0^1 g(y) dy. y(x) = sigma(2c), c = (pi/2) sinh(x),
# sigma = logistic. y, log(y) and -log(1-y) via stats::plogis(..., log.p=) stay stable for large
# |c| (no underflow of y to exactly 0/1, which would make y^m = Inf at far nodes).

#' Tanh-Sinh Quadrature Rule on (0, 1)
#'
#' Computes tanh-sinh (double-exponential) quadrature nodes and weights for
#' integrals over \eqn{[0, 1]}. The transformation \eqn{y = \sigma(2c)} with
#' \eqn{c = (\pi/2) \sinh(x)} maps the real line to \eqn{(0, 1)} with doubly
#' exponential endpoint decay. Log-scale values \code{logy} and \code{-log(1-y)}
#' are returned for stable evaluation of endpoint-singular integrands.
#'
#' @param h Node spacing on the transformed axis. Default \code{1/32}.
#' @param X Half-width on the transformed axis. Default \code{3}.
#'
#' @return A list with elements \code{y} (nodes in \eqn{(0, 1)}), \code{w}
#'   (weights), \code{logy} (\eqn{\log(y)}), and \code{g0}
#'   (\eqn{-\log(1-y)}).
#' @noRd
.tanh_sinh <- function(h = 1/32, X = 3) {
  x  <- seq(-X, X, by = h)
  c. <- (pi / 2) * sinh(x)
  y    <- stats::plogis(2 * c.)
  logy <- stats::plogis(2 * c., log.p = TRUE)                       # log(y), stable
  g0   <- -stats::plogis(2 * c., lower.tail = FALSE, log.p = TRUE)  # -log(1-y), stable
  w  <- h * (pi / 2) * cosh(x) / (2 * cosh(c.)^2)
  list(y = y, w = w, logy = logy, g0 = g0)
}

#' Cache Environment for Tanh-Sinh Rules
#'
#' Environment storing pre-computed tanh-sinh quadrature rules keyed by
#' \code{(h, X)} pair. Populated lazily by \code{.ts_rule}.
#'
#' @noRd
.ts_cache <- new.env(parent = emptyenv())

#' Cached Tanh-Sinh Quadrature Rule
#'
#' Returns a tanh-sinh quadrature rule for the given node spacing and
#' half-width, retrieving from cache if available or computing via
#' \code{.tanh_sinh} and caching the result.
#'
#' @param h Node spacing on the transformed axis. Default \code{1/32}.
#' @param X Half-width on the transformed axis. Default \code{3}.
#'
#' @return A list with elements \code{y}, \code{w}, \code{logy}, \code{g0}.
#' @noRd
.ts_rule <- function(h = 1/32, X = 3) {
  key <- sprintf("%g_%g", h, X)
  if (is.null(.ts_cache[[key]])) .ts_cache[[key]] <- .tanh_sinh(h, X)
  .ts_cache[[key]]
}

#' Quadrature L-Moments of the Exponentiated Weibull Distribution
#'
#' Computes L-moments of the exponentiated Weibull distribution via tanh-sinh
#' (double-exponential) quadrature of the probability-weighted moments. The
#' tanh-sinh rule ensures doubly exponential convergence at both endpoints,
#' giving machine-precision accuracy for all positive parameter values with
#' pure \eqn{\mathcal{O}(N)} arithmetic — no adaptive integration required.
#' All L-moments exist for any \eqn{a, b, s > 0}.
#'
#' @param order Integer vector of L-moment orders (e.g. \code{1:5}).
#' @param scale Scale parameter (\eqn{s > 0}).
#' @param shape1 First shape parameter (\eqn{a > 0}), the Weibull shape;
#'   \eqn{a = 1} gives the exponential baseline.
#' @param shape2 Second shape parameter (\eqn{b > 0}), the exponentiation
#'   power; \eqn{b = 1} recovers the ordinary Weibull.
#' @param h Tanh-sinh step size. Automatically chosen if \code{NULL}.
#' @param X Tanh-sinh half-width. Automatically chosen if \code{NULL}.
#'
#' @return Named numeric vector of L-moments (\code{lambda_1}, \code{lambda_2},
#'   \code{tau_3}, \code{tau_4}, ...).
#' @keywords internal
#' @noRd
lmom_expweibull <- function(order = 1:5, scale, shape1, shape2, h = NULL, X = NULL) {
  a <- shape1; b <- shape2; s <- scale
  rmax <- max(order)
  # Endpoint singularity strength at y=0 (integrand ~ y^(b-1+1/a) at r=0): alpha < 0 (large shape1
  # & small shape2) needs nodes reaching closer to 0; alpha >= 0 is benign (all realistic rainfall).
  if (is.null(h) || is.null(X)) {
    if (b - 1 + 1 / a < 0) { h <- 1/64; X <- 5 } else { h <- 1/32; X <- 3.5 }
  }
  ts <- .ts_rule(h, X)
  base <- ts$w * ts$g0^(1 / a)            # w * (-log(1-y))^(1/a)  (depends on a only)

  beta <- vapply(0:(rmax - 1), function(r)
    s * b * sum(base * exp((b * (r + 1) - 1) * ts$logy)), numeric(1))

  lambda <- .pwm_to_lmom(beta)
  stats::setNames(ifelse(order <= 2, lambda[order], lambda[order] / lambda[2]),
                  ifelse(order <= 2, paste0("lambda_", order), paste0("tau_", order)))
}
