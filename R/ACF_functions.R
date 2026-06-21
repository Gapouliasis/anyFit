#' Theoretical autocorrelation functions
#'
#' @description
#' Closed-form expressions for three commonly used theoretical autocorrelation function (ACF)
#' models: the Cauchy-type algebraic decay
#' (CAS), the Hurst--Kolmogorov ACF (HK), and the
#' short-range dependence exponential decay (SRD). Each function returns a
#' data frame of ACF values for lags 0 through \code{lag_max}.
#'
#' @name ACF_functions
#' @export

#' @rdname ACF_functions
#'
#' @description
#' Cauchy-type autocorrelation structure with algebraic decay, governed by two
#' parameters: \code{kappa} controls the scale and \code{beta} controls the
#' decay rate. The ACF decays hyperbolically, producing long-range dependence
#' for small \eqn{\beta}.
#'
#' The ACF is:
#' \deqn{\rho(\tau) = (1 + \kappa\beta\tau)^{-1/\beta}}{\rho(\tau) = (1 + k*b*t)^(-1/b)}
#' where:
#' * \eqn{\kappa}{k} --- scale parameter (\code{kappa})
#' * \eqn{\beta}{b} --- shape parameter (\code{beta}), controls the algebraic decay rate
#' * \eqn{\tau}{t} --- lag
#'
#' @param kappa Scale parameter. Positive numeric scalar.
#' @param beta Shape parameter. Controls the algebraic decay rate. Positive numeric scalar.
#' @param lag_max Maximum lag for which ACF values are computed. Default is 10.
#'
#' @return A data frame with columns \code{lag} (integer, 0 to \code{lag_max})
#'   and \code{ACF} (numeric, the theoretical autocorrelation at each lag).
#'
#' @examples
#' # Cauchy-type ACF with slow decay
#' CAS_ACF(kappa = 0.5, beta = 0.3)
#'
#' # Faster decay with larger beta
#' CAS_ACF(kappa = 0.5, beta = 1.5, lag_max = 20)
#'
#' @export
CAS_ACF <- function(kappa,beta,lag_max = 10){
  p_cas <- data.frame(lag = seq(0,lag_max))
  p_cas$ACF <- (1+kappa*beta*p_cas$lag)^(-1/beta)
  return(p_cas)
}

#' @rdname ACF_functions
#'
#' @description
#' Hurst--Kolmogorov autocorrelation function, governed by a single Hurst exponent \eqn{H}. Values
#' \eqn{H > 0.5} indicate long-range dependence, \eqn{H = 0.5} corresponds to
#' white noise, and \eqn{H < 0.5} indicates anti-persistence.
#'
#' The ACF is:
#' \deqn{\rho(\tau) = \frac{1}{2}\left(|\tau-1|^{2H} - 2|\tau|^{2H} + |\tau+1|^{2H}\right)}{\rho(\tau) = 0.5 * (|t-1|^(2H) - 2*|t|^(2H) + |t+1|^(2H))}
#' where:
#' * \eqn{H}{H} --- Hurst exponent (\code{H}), \eqn{0 < H < 1}{0 < H < 1}; \eqn{H > 0.5}{H > 0.5} indicates long-range dependence, \eqn{H = 0.5}{H = 0.5} is white noise, \eqn{H < 0.5}{H < 0.5} indicates anti-persistence
#' * \eqn{\tau}{t} --- lag
#'
#' @param H Hurst exponent. Numeric scalar in \eqn{(0, 1)}.
#' @param lag_max Maximum lag for which ACF values are computed. Default is 10.
#'
#' @return A data frame with columns \code{lag} (integer, 0 to \code{lag_max})
#'   and \code{ACF} (numeric, the theoretical autocorrelation at each lag).
#'
#' @examples
#' # Long-range dependence (H > 0.5)
#' HK_ACF(H = 0.8)
#'
#' # Anti-persistence (H < 0.5)
#' HK_ACF(H = 0.3, lag_max = 20)
#'
#' @export
HK_ACF <- function(H , lag_max = 10){
  p_hk = data.frame(lag = seq(0,lag_max))
  p_hk$ACF <- 0.5*((abs(p_hk$lag-1))^(2*H)-
                     2*(abs(p_hk$lag))^(2*H)+
                     (abs(p_hk$lag+1))^(2*H))
  return(p_hk)
}

#' @rdname ACF_functions
#'
#' @description
#' Short-range dependence autocorrelation function with exponential decay,
#' corresponding to a first-order Markov (AR1) process. Governed by a single
#' decay parameter \eqn{\kappa}, the lag-1 autocorrelation is
#' \eqn{\rho_1 = \exp(-\kappa)}.
#'
#' The ACF is:
#' \deqn{\rho(\tau) = \exp(-\kappa\tau)}{\rho(\tau) = exp(-k*t)}
#' where:
#' * \eqn{\kappa}{k} --- decay rate (\code{kappa}); the lag-1 autocorrelation is \eqn{\rho_1 = \exp(-\kappa)}{\rho_1 = exp(-k)}
#' * \eqn{\tau}{t} --- lag
#'
#' @param kappa Decay rate. Positive numeric scalar.
#' @param lag_max Maximum lag for which ACF values are computed. Default is 10.
#'
#' @return A data frame with columns \code{lag} (integer, 0 to \code{lag_max})
#'   and \code{ACF} (numeric, the theoretical autocorrelation at each lag).
#'
#' @examples
#' # Moderate short-range dependence
#' SRD_ACF(kappa = 0.5)
#'
#' # Very short memory
#' SRD_ACF(kappa = 2, lag_max = 15)
#'
#' @export
SRD_ACF <- function(kappa, lag_max = 10){
  psrd = data.frame(lag = seq(0,lag_max))
  psrd$ACF <- exp(-psrd$lag*kappa)
  return(psrd)
}
