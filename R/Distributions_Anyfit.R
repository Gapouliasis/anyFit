# Two-sample Cramer-von Mises statistic computed on already ascending-sorted
# inputs. Numerically identical to CDFt::CramerVonMisesTwoSamples
# (verified including ties).
cvm_two_sample_sorted = function(xS1, xS2){
  M <- length(xS1); N <- length(xS2)
  somM <- sum(findInterval(xS1, xS2, left.open = TRUE)^2)  # #{xS2 <  xS1}
  somN <- sum(findInterval(xS2, xS1)^2)                    # #{xS1 <= xS2}
  U <- N * somN + M * somM
  ((U / (N * M)) / (N + M)) - ((4 * M * N - 1) / (6 * (M + N)))
}

# Two-sample Cramer-von Mises statistic. Wrapper that sorts inputs and
# delegates to the pre-sorted engine.
cvm_two_sample = function(S1, S2) cvm_two_sample_sorted(sort(S1), sort(S2))

# Two-sample Kolmogorov-Smirnov statistic computed on already ascending-sorted
# inputs. Identical to CDFt::KolmogorovSmirnov.
ks_two_sample_sorted = function(xS1, xS2){
  M <- length(xS1); N <- length(xS2)
  cdf1     <- findInterval(xS1, xS1) / M
  cdfEstim <- findInterval(xS2, xS2) / N
  cdfRef   <- stats::approx(xS1, cdf1, xS2, yleft = 0, yright = 1, ties = "mean")$y
  max(abs(cdfRef - cdfEstim))
}

# Two-sample Kolmogorov-Smirnov statistic. Wrapper that sorts inputs and
# delegates to the pre-sorted engine.
ks_two_sample = function(S1, S2) ks_two_sample_sorted(sort(S1), sort(S2))

# Computes six goodness-of-fit metrics for a fitted distribution against
# empirical data: maximum-likelihood criterion (MLE), two-sample
# Cramer-von Mises (CM), two-sample Kolmogorov-Smirnov (KS), mean
# squared quantile error (MSEquant), percentage difference of maxima
# (DiffOfMax), and mean absolute difference of the ten largest values
# (MeanDiffOf10Max). Both the empirical and fitted quantile vectors are
# sorted once and reused across all metrics. The function matches the
# d/p/q/r family by prefix (e.g. distribution = "exp" resolves to dexp,
# pexp, qexp, rexp).
GOF_tests = function(x, fit, distribution){
  qfun <- match.fun(paste0('q', distribution))
  dfun <- match.fun(paste0('d', distribution))
  q_emp <- coredata(x)
  n <- length(q_emp)
  u_emp <- rank(q_emp, na.last = NA, ties.method = "average")/(n+1)
  qq_fitted <- do.call(qfun, c(list(p = u_emp), fit))

  # Sort each vector once and reuse for CM, KS, the top-10 block and the max.
  xq <- sort(q_emp)
  xf <- sort(qq_fitted)

  CM <- cvm_two_sample_sorted(xq, xf)
  KS <- ks_two_sample_sorted(xq, xf)

  MLE = -sum(log(do.call(dfun, c(list(x = x), fit)) + 0.000001))

  MSEquant<-sum((qq_fitted-q_emp)^2)/n
  DiffOfMax<-((xf[n]-xq[n])/xq[n])*100
  MeanDiffOf10Max<-sum(abs(rev(xf)[1:10]-rev(xq)[1:10]))/10

  # AIC = 2*length(fit)-2*sum(log(do.call(paste0('d',distribution),c(list(x = x), fit))))
  # BIC = length(fit)*log(length(x))-2*sum(do.call(paste0('d',distribution),c(list(x = x), fit)))

  GoF <- list(MLE=MLE, CM = CM, KS = KS,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  return(GoF)
}


# Negative log-likelihood objective for maximum-likelihood parameter
# optimisation. Accepts a trial parameter vector, an xts data series, and
# the density function of a distribution, and returns the negative summed
# log-likelihood. Infinite values are replaced with a large finite penalty
# (1e10) to keep gradient-based optimisers from diverging.
MLE_fun = function(trial = c(par1,par2,par3),x_ts,dfunction){
  #trial = c(par1,par2,par3)
  dist_args = formalArgs(dfunction)
  temp_args = list(x_ts)
  temp_args = c(temp_args, as.list(trial),1)
  names(temp_args) = dist_args
  LF = do.call(dfunction, temp_args)
  LF = -sum(log(LF))
  if (is.infinite(LF)){LF = 10^10}
  return(LF)
}


# Exponential -------------------------------------------------------------

#' @name Exponential
#' @title Exponential Distribution
#'
#' @description Two-parameter Exponential distribution with location
#'   \eqn{\mu}{m} and scale \eqn{\sigma}{s}. The distribution is
#'   lower-bounded at the location parameter. Fitting
#'   is performed via closed-form L-moments. An optional fixed lower bound
#'   \code{bound} is supported for cases where the minimum is known a
#'   priori.
#'
#'   The probability density function is:
#'   \deqn{f(x) = \frac{1}{\sigma} \exp\left(-\frac{x-\mu}{\sigma}\right), \quad x \ge \mu}{f(x) = (1/s) exp(-(x-m)/s), x >= m}
#'   where:
#'   * \eqn{\mu}{m} --- location parameter (\code{location})
#'   * \eqn{\sigma}{s} --- scale parameter (\code{scale})
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @export
#'

# x<-rexp(1000,location=0,scale=3) # my function

dexp=function(x,location,scale) {
  fx=stats::dexp(x=x-location,rate=1/scale)
  return(fx)
}

#' @rdname Exponential
#' @export
pexp=function(q,location,scale) {
  FX=lmom::cdfexp(x=q,para=c(location,scale))
  return(FX)
}
#' @rdname Exponential
#' @export
qexp=function(p,location,scale) {
  x=lmom::quaexp(f=p,para=c(location,scale))
  return(x)
}
#' @rdname Exponential
#' @export
rexp=function(n,location,scale) {
  x=lmom::quaexp(f=stats::runif(n),para=c(location,scale))
  return(x)
}


#' @title fitlm_exp
#'
#' @description Fits the two-parameter Exponential distribution to an xts
#'   series by the method of L-moments. If \code{bound} is supplied, the
#'   location is fixed to that value and the scale is derived from the
#'   second L-moment; otherwise both parameters are obtained from
#'   \code{\link[lmom]{pelexp}}. Zero values below \code{zero_threshold}
#'   may be excluded via \code{ignore_zeros}. Goodness-of-fit is assessed
#'   via \code{GOF_tests}, and theoretical L-moments are computed from the
#'   fitted parameters for comparison with the sample.
#'
#' @param x An xts object containing the time series data.
#' @param bound Numeric or \code{NULL}. Optional fixed lower bound
#'   (location). Default \code{NULL}.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted parameters), \code{TheorLMom} (theoretical L-moments),
#'   \code{DataLMom} (sample L-moments), and \code{GoF} (goodness-of-fit
#'   metrics).
#'
#' @examples
#' x <- xts::xts(rexp(365, location = 0, scale = 10), order.by = as.Date("2020-01-01") + 0:364)
#' fit <- fitlm_exp(x)
#' fit$Param
#'
#' @export
#'

fitlm_exp=function(x,bound=NULL, ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- stats::na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }
  sam=lmom::samlmu(x)

  if(is.null(bound)){

    par=lmom::pelexp(sam)
    names(par) = NULL
    fit = list()
    fit$location=par[1]
    fit$scale=par[2]

  } else {

    names(sam)=NULL
    fit=list()
    fit$location=bound
    fit$scale=sam[2]*2
    par<-c(fit$location,fit$scale)

  }

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'exp')

  Res<-list()
  Res$Distribution<-list(FXs="qexp",bound=bound)
  Res$Param<-fit
  Res$TheorLMom<-lmom::lmrexp(par,nmom=4)
  names(sam)<-c('lambda_1','lambda_2','tau_3','tau_4')
  Res$DataLMom<-sam
  Res$GoF<-GoF

  return(Res)
}


# Rayleigh ----------------------------------------------------------
#' @name Rayleigh
#' @title Rayleigh Distribution
#'
#' @description Two-parameter Rayleigh distribution with location
#'   \eqn{\mu}{m} and scale \eqn{\sigma}{s}. The distribution is
#'   lower-bounded at the location parameter. Fitting uses closed-form L-moment expressions
#'   derived from gamma-function relationships; the location is obtained
#'   from the first two L-moments and the scale from the second L-moment
#'   ratio.
#'
#'   The probability density function is:
#'   \deqn{f(x) = \frac{x - \mu}{\sigma^2} \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right), \quad x \ge \mu}{f(x) = (x-m)/s^2 * exp(-(x-m)^2 / (2s^2)), x >= m}
#'   where:
#'   * \eqn{\mu}{m} --- location parameter (\code{location})
#'   * \eqn{\sigma}{s} --- scale parameter (\code{scale})
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @export
#'
drayleigh = function(x, location, scale) {
  fx <- VGAM::drayleigh(x - location, scale, log = FALSE)
  return(fx)
}

#' @rdname Rayleigh
#' @export
prayleigh = function(x, location, scale) {
  fx <- VGAM::prayleigh(x - location, scale, lower.tail = TRUE, log.p = FALSE)
  return(fx)
}
#' @rdname Rayleigh
#' @export
qrayleigh = function(p, location, scale) {
  fx <- VGAM::qrayleigh(p, scale, lower.tail = TRUE, log.p = FALSE) + location
  return(fx)
}
#' @rdname Rayleigh
#' @export
rrayleigh = function(n, location, scale) {
  fx <- VGAM::rrayleigh(n, scale) + location
  return(fx)
}


#' @title fitlm_rayleigh
#'
#' @description Fits the two-parameter Rayleigh distribution to an xts
#'   series by the method of L-moments. Closed-form expressions derived
#'   from gamma-function integrals are used: the location is recovered
#'   from \eqn{\lambda_1} and \eqn{\lambda_2}, and the scale follows from
#'   the squared second L-moment ratio. Zero values below
#'   \code{zero_threshold} may be excluded via \code{ignore_zeros}.
#'   Goodness-of-fit is assessed via \code{GOF_tests}.
#'
#' @param x An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted parameters), \code{TheorLMom} (theoretical L-moments),
#'   \code{DataLMom} (sample L-moments), and \code{GoF} (goodness-of-fit
#'   metrics).
#'
#' @examples
#' x <- xts::xts(rrayleigh(365, location = 0, scale = 2),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' fit <- fitlm_rayleigh(x)
#' fit$Param
#'
#' @export
#'

fitlm_rayleigh = function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- stats::na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }

  sam=lmom::samlmu(x)
  l1 <- lmom::samlmu(x,nmom = 1, ratios = FALSE)
  l2 <- lmom::samlmu(x,nmom = 2, ratios = FALSE)
  location <- l1 - sqrt(2)/(sqrt(2) - 1)*l2
  ratio <- (gamma(1.5))^2*(3-2*sqrt(2))/(2*l2^2)
  scale <- 1/sqrt(2*ratio)

  fit = list(location = location, scale = scale)

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'rayleigh')

  Res<-list()
  Res$Distribution<-list(FXs="qrayleigh")
  Res$Param<-fit
  Res$TheorLMom<-c(gamma(1.5)/sqrt(scale)+location, gamma(1.5)/sqrt(scale)*(sqrt(2)-1)/sqrt(2))
  names(sam)<-c('lambda_1','lambda_2','tau_3','tau_4')
  Res$DataLMom<-sam
  Res$GoF<-GoF

  return(Res)
}



# Gamma -------------------------------------------------------------

# x<-stats::rgamma(10000,scale=3,shape=0.5)

#' @name Gamma
#' @title Gamma Distribution
#'
#' @description Two-parameter Gamma distribution with scale \eqn{\beta}{b}
#'   and shape \eqn{\alpha}{a}. The distribution is supported on
#'   \eqn{x > 0} and is widely used for modelling positively skewed
#'   variables such as precipitation amounts. Density, distribution,
#'   quantile, and random generation functions are provided via base R
#'   \code{\link[stats]{dgamma}}. Fitting is performed via closed-form
#'   L-moments through \code{\link[lmom]{pelgam}}.
#'
#'   The probability density function is:
#'   \deqn{f(x) = \frac{1}{\beta^\alpha \Gamma(\alpha)} x^{\alpha-1} \exp\left(-\frac{x}{\beta}\right), \quad x > 0}{f(x) = x^(a-1) * exp(-x/b) / (b^a * Gamma(a)), x > 0}
#'   where:
#'   * \eqn{\beta}{b} --- scale parameter (\code{scale})
#'   * \eqn{\alpha}{a} --- shape parameter (\code{shape})
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param scale scale parameter.
#' @param shape shape parameter.
#' @export
#'

dgamma=function(x,scale,shape) {
  fx=stats::dgamma(x=x,shape=shape,scale=scale)
  return(fx)
}

#' @rdname Gamma
#' @export
pgamma=function(q,scale,shape) {
  FX=stats::pgamma(q=q,shape=shape,scale=scale)
  return(FX)
}

#' @rdname Gamma
#' @export
qgamma=function(p,scale,shape) {
  x=stats::qgamma(p=p,shape=shape,scale=scale)
  return(x)
}

#' @rdname Gamma
#' @export
rgamma=function(n,scale,shape) {
  x=stats::rgamma(n=n,shape=shape,scale=scale)
  return(x)
}

#' @title fitlm_gamma
#'
#' @description Fits the two-parameter Gamma distribution to an xts series
#'   by the method of L-moments. Parameters are obtained in closed form
#'   via \code{\link[lmom]{pelgam}}, which matches the sample L-moments
#'   to the theoretical L-moment ratios of the Gamma distribution. Zero
#'   values below \code{zero_threshold} may be excluded via
#'   \code{ignore_zeros}. Goodness-of-fit is assessed via
#'   \code{GOF_tests}, and theoretical L-moments are computed from the
#'   fitted parameters for comparison with the sample.
#'
#' @param x An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted \code{scale} and \code{shape}), \code{TheorLMom}
#'   (theoretical L-moments), \code{DataLMom} (sample L-moments), and
#'   \code{GoF} (goodness-of-fit metrics).
#'
#' @examples
#' x <- xts::xts(rgamma(365, scale = 3, shape = 0.5),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' fit <- fitlm_gamma(x)
#' fit$Param
#'
#' @export
#'

fitlm_gamma=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- stats::na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }

  sam=lmom::samlmu(x)
  par=lmom::pelgam(sam)
  names(par)=NULL
  fit=list()
  fit$scale = par[2]
  fit$shape = par[1]

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'gamma')

  Res<-list()
  Res$Distribution<-list(FXs="qgamma")
  Res$Param<-fit
  Res$TheorLMom<-lmom::lmrgam(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF<-GoF


  return(Res)
}


# Gamma3 -------------------------------------------------------------

# x<-rgamma3(10000,location=5,scale=0.5,shape=3)
#' @name Gamma3
#' @title 3-parameter Gamma Distribution
#'
#' @description Three-parameter Gamma (Pearson type III) distribution with
#'   location \eqn{\mu}{m}, scale \eqn{\beta}{b}, and shape
#'   \eqn{\alpha}{a}. The distribution is lower-bounded at the location
#'   parameter and is widely employed in flood-frequency analysis. Density,
#'   distribution, quantile, and random generation functions are provided
#'   via the \pkg{PearsonDS} package. Fitting is performed via closed-form
#'   L-moments through \code{\link[lmom]{pelpe3}}.
#'
#'   The probability density function is:
#'   \deqn{f(x) = \frac{\beta^\alpha}{\Gamma(\alpha)} (x - \mu)^{\alpha-1} \exp\left(-\beta (x - \mu)\right), \quad x \ge \mu}{f(x) = b^a / Gamma(a) * (x-m)^(a-1) * exp(-b*(x-m)), x >= m}
#'   where:
#'   * \eqn{\mu}{m} --- location parameter (\code{location})
#'   * \eqn{\beta}{b} --- scale parameter (\code{scale})
#'   * \eqn{\alpha}{a} --- shape parameter (\code{shape})
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param shape shape parameter.
#' @export
#'

dgamma3=function(x,location,scale,shape) {
  fx=PearsonDS::dpearsonIII(x=x,shape=shape,location=location,scale=1/scale)
  return(fx)
}

#' @rdname Gamma3
#' @export
pgamma3=function(q,location,scale,shape) {
  FX=PearsonDS::ppearsonIII(q=q,shape=shape,location=location,scale=1/scale)
  return(FX)
}

#' @rdname Gamma3
#' @export
qgamma3=function(p,location,scale,shape) {
  x=PearsonDS::qpearsonIII(p=p,shape=shape,location=location,scale=1/scale)
  return(x)
}

#' @rdname Gamma3
#' @export
rgamma3=function(n,location,scale,shape) {
  x=PearsonDS::rpearsonIII(n=n,shape=shape,location=location,scale=1/scale)
  return(x)
}

#' @title fitlm_gamma3
#'
#' @description Fits the three-parameter Gamma (Pearson type III)
#'   distribution to an xts series by the method of L-moments. Parameters
#'   are obtained via \code{\link[lmom]{pelpe3}} and then transformed to
#'   the standard location, scale, and shape parameterisation used
#'   throughout the package. Zero values below \code{zero_threshold} may be
#'   excluded via \code{ignore_zeros}. Goodness-of-fit is assessed via
#'   \code{GOF_tests}, and theoretical L-moments are computed from the
#'   fitted parameters.
#'
#' @param x An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted parameters), \code{TheorLMom} (theoretical L-moments),
#'   \code{DataLMom} (sample L-moments), and \code{GoF} (goodness-of-fit
#'   metrics).
#'
#' @examples
#' x <- xts::xts(rgamma3(365, location = 5, scale = 0.5, shape = 3),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' fit <- fitlm_gamma3(x)
#' fit$Param
#'
#' @export
#'


fitlm_gamma3=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- stats::na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }

  sam=lmom::samlmu(x)
  par=lmom::pelpe3(sam)
  names(par)=NULL

  fit=list()
  fit$location=par[1]-2*par[2]/par[3]
  fit$scale=1/((1/2)*par[2]*abs(par[3]))
  fit$shape=4/(par[3]^2)

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'gamma3')

  Res<-list()
  Res$Distribution<-list(FXs="qgamma3")
  Res$Param<-fit
  Res$TheorLMom<-lmom::lmrpe3(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF<-GoF

  return(Res)

}


# Generalised Logistic -------------------------------------------------------------

# x<-quaglo(runif(10000), c(0,1,-0.5))

#' @name GenLogistic
#' @title Generalized Logistic Distribution
#'
#' @description Three-parameter Generalised Logistic distribution with
#'   location \eqn{\mu}{m}, scale \eqn{\sigma}{s}, and shape
#'   \eqn{\xi}{j}. Density, distribution, and quantile functions are
#'   provided via the \pkg{lmomco} package. The distribution encompasses
#'   symmetric (\eqn{\xi = 0}) and asymmetric forms and is used in
#'   regional frequency analysis of hydrological extremes. Fitting is
#'   performed via closed-form L-moments through
#'   \code{\link[lmom]{pelglo}}.
#'
#'   Define \eqn{z = (x - \mu)/\sigma}{z = (x-m)/s}. The probability
#'   density function is:
#'   \deqn{f(x) = \frac{1}{\sigma} \frac{e^{-z}}{(1+e^{-z})^2}, \quad \xi = 0}{f(x) = exp(-z) / (s*(1+exp(-z))^2)}
#'   \deqn{f(x) = \frac{1}{\sigma} \frac{(1-\xi z)^{1/\xi-1}}{(1+(1-\xi z)^{1/\xi})^2}, \quad \xi \neq 0}{f(x) = (1-j*z)^(1/j-1) / (s*(1+(1-j*z)^(1/j))^2)}
#'   where:
#'   * \eqn{\mu}{m} --- location parameter (\code{location})
#'   * \eqn{\sigma}{s} --- scale parameter (\code{scale})
#'   * \eqn{\xi}{j} --- shape parameter (\code{shape})
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param shape shape parameter.
#' @export
#'
dgenlogi <- function(x, location, scale, shape) {
  para <- lmomco::vec2par(c(location, scale, shape), type = "glo")
  lmomco::pdfglo(x, para)
}
#' @rdname GenLogistic
#' @export
pgenlogi <- function(q, location, scale, shape) {
  para <- lmomco::vec2par(c(location, scale, shape), type = "glo")
  lmomco::cdfglo(q, para)
}
#' @rdname GenLogistic
#' @export
qgenlogi <- function(p, location, scale, shape) {
  para <- lmomco::vec2par(c(location, scale, shape), type = "glo")
  lmomco::quaglo(p, para)
}

#' @title fitlm_genlogi
#'
#' @description Fits the three-parameter Generalised Logistic distribution
#'   to an xts series by the method of L-moments. Parameters are obtained
#'   in closed form via \code{\link[lmom]{pelglo}}, which matches the
#'   sample L-moment ratios to the theoretical L-moment ratios of the
#'   Generalised Logistic distribution. Zero values below
#'   \code{zero_threshold} may be excluded via \code{ignore_zeros}.
#'   Goodness-of-fit is assessed via \code{GOF_tests}, and theoretical
#'   L-moments are computed from the fitted parameters for comparison with
#'   the sample.
#'
#' @param x An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted parameters), \code{TheorLMom} (theoretical L-moments),
#'   \code{DataLMom} (sample L-moments), and \code{GoF} (goodness-of-fit
#'   metrics).
#'
#' @examples
#' x <- xts::xts(qgenlogi(runif(365), location = 0, scale = 1, shape = -0.5),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' fit <- fitlm_genlogi(x)
#' fit$Param
#'
#' @export
#'


fitlm_genlogi=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- stats::na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }

  sam=lmom::samlmu(x)
  par=lmom::pelglo(sam)
  names(par)=NULL
  fit=list()
  fit$location = par[1]
  fit$scale = par[2]
  fit$shape = par[3]

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'genlogi')

  Res<-list()
  Res$Distribution<-list(FXs="qgenlogi")
  Res$Param<-fit
  Res$TheorLMom<-lmom::lmrglo(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF<-GoF

  return(Res)
}


# Normal -------------------------------------------------------------

# x<-quanor(runif(10000), c(0,3))
# x<-rnorm(10000, mean=0,sd=3)
# identical runs

#' @name Normal
#' @title Normal Distribution
#'
#' @description Two-parameter Normal (Gaussian) distribution with mean
#'   \eqn{\mu}{m} and standard deviation \eqn{\sigma}{s}. Fitting is performed via
#'   closed-form L-moments. The
#'   L-moment estimates coincide with the conventional method-of-moments
#'   estimates (\code{mean} = \eqn{\lambda_1}, \code{sd} =
#'   \eqn{\lambda_2 \sqrt{\pi}}).
#'
#'   The probability density function is:
#'   \deqn{f(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{1}{2}\left(\frac{x - \mu}{\sigma}\right)^2\right)}{f(x) = exp(-0.5*((x-m)/s)^2) / (s*sqrt(2*pi))}
#'   where:
#'   * \eqn{\mu}{m} --- location/mean parameter (\code{mean})
#'   * \eqn{\sigma}{s} --- scale/standard deviation parameter (\code{sd})
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param mean mean parameter.
#' @param sd standard deviation parameter.
#' @export
#'

dnorm=function(x,mean=0,sd=1) {
  fx=stats::dnorm(x=x,sd=sd,mean=mean)
  return(fx)
}

#' @rdname Normal
#' @export
pnorm=function(q,mean=0,sd=1) {
  FX=stats::pnorm(q=q,sd=sd,mean=mean)
  return(FX)
}

#' @rdname Normal
#' @export
qnorm=function(p,mean=0,sd=1) {
  x=stats::qnorm(p=p,sd=sd,mean=mean)
  return(x)
}

#' @rdname Normal
#' @export
rnorm=function(n,mean=0,sd=1) {
  x=stats::rnorm(n=n,sd=sd,mean=mean)
  return(x)
}

#' @title fitlm_norm
#'
#' @description Fits the Normal distribution to an xts series by the
#'   method of L-moments. Parameters are obtained in closed form via
#'   \code{\link[lmom]{pelnor}}, which maps the first two sample L-moments
#'   directly to the mean and standard deviation. Zero values below
#'   \code{zero_threshold} may be excluded via \code{ignore_zeros}.
#'   Goodness-of-fit is assessed via \code{GOF_tests}, and theoretical
#'   L-moments are computed from the fitted parameters for comparison with
#'   the sample.
#'
#' @param x An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted \code{mean} and \code{sd}), \code{TheorLMom}
#'   (theoretical L-moments), \code{DataLMom} (sample L-moments), and
#'   \code{GoF} (goodness-of-fit metrics).
#'
#' @examples
#' x <- xts::xts(rnorm(365, mean = 0, sd = 3),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' fit <- fitlm_norm(x)
#' fit$Param
#'
#' @export
#'


fitlm_norm=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- stats::na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }

  sam=lmom::samlmu(x)
  par=lmom::pelnor(sam)
  names(par) = NULL
  fit = list()
  fit$mean=par[1]
  fit$sd = par[2]

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'norm')

  Res<-list()
  Res$Distribution<-list(FXs="qnorm")
  Res$Param<-fit
  Res$TheorLMom<-lmom::lmrnor(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF<-GoF

  return(Res)
}



# Weibull, 3-parameter -------------------------------------------------------------

# x<-quawei(runif(10000), c(0,5,1.5)) # require(lmom)
# x<-rweibull(10,scale = 5,shape = 1.5)
# identical runs

#' @name Weibull3
#' @title 3-parameter Weibull Distribution
#'
#' @description Three-parameter Weibull distribution with location
#'   \eqn{\mu}{m}, scale \eqn{\sigma}{s}, and shape \eqn{k}{k}. The
#'   distribution is lower-bounded at the location parameter and
#'   generalises the Exponential (\eqn{k = 1}) and Rayleigh (\eqn{k = 2})
#'   distributions. Fitting is
#'   performed via closed-form L-moments. An optional fixed lower bound is
#'   supported.
#'
#'   The probability density function is:
#'   \deqn{f(x) = \frac{k}{\sigma} \left(\frac{x - \mu}{\sigma}\right)^{k-1} \exp\left(-\left(\frac{x - \mu}{\sigma}\right)^k\right), \quad x \ge \mu}{f(x) = k/s * ((x-m)/s)^(k-1) * exp(-((x-m)/s)^k), x >= m}
#'   where:
#'   * \eqn{\mu}{m} --- location parameter (\code{location})
#'   * \eqn{\sigma}{s} --- scale parameter (\code{scale})
#'   * \eqn{k}{k} --- shape parameter (\code{shape})
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param shape shape parameter.
#' @export
#'

dweibull=function(x,location,scale,shape) {
  fx=FAdist::dweibull3(x=x,shape=shape,scale=scale,thres=location)
  return(fx)
}

#' @rdname Weibull3
#' @export
pweibull=function(q,location,scale,shape) {
  FX=FAdist::pweibull3(q=q,shape=shape,scale=scale,thres=location)
  return(FX)
}

#' @rdname Weibull3
#' @export
qweibull=function(p,location,scale,shape) {
  x=FAdist::qweibull3(p=p,shape=shape,scale=scale,thres=location)
  return(x)
}

#' @rdname Weibull3
#' @export
rweibull=function(n,location,scale,shape) {
  x=FAdist::rweibull3(n=n,shape=shape,scale=scale,thres=location)
  return(x)
}

#' @title fitlm_weibull
#'
#' @description Fits the three-parameter Weibull distribution to an xts
#'   series by the method of L-moments. Parameters are obtained in closed
#'   form via \code{\link[lmom]{pelwei}}, which matches the sample
#'   L-moment ratios to the theoretical L-moment ratios of the Weibull
#'   distribution. If \code{bound} is supplied, the location is fixed to
#'   that value; otherwise all three parameters are free. Zero values
#'   below \code{zero_threshold} may be excluded via
#'   \code{ignore_zeros}. Goodness-of-fit is assessed via
#'   \code{GOF_tests}.
#'
#' @param x An xts object containing the time series data.
#' @param bound Numeric or \code{NULL}. Optional fixed lower bound
#'   (location). Default \code{NULL}.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted parameters), \code{TheorLMom} (theoretical L-moments),
#'   \code{DataLMom} (sample L-moments), and \code{GoF} (goodness-of-fit
#'   metrics).
#'
#' @examples
#' x <- xts::xts(rweibull(365, location = 0, scale = 5, shape = 1.5),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' fit <- fitlm_weibull(x)
#' fit$Param
#'
#' @export
#'

fitlm_weibull=function(x,bound=NULL,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- stats::na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }

  sam=lmom::samlmu(x)
  par=lmom::pelwei(sam,bound=bound)
  names(par) = NULL
  fit=list()
  fit$location=par[1]
  fit$scale = par[2]
  fit$shape = par[3]

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'weibull')

  Res<-list()
  Res$Distribution<-list(FXs="qweibull",bound=bound)
  Res$Param<-fit
  Res$TheorLMom<-lmom::lmrwei(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF<-GoF

  return(Res)
}


# Gumbel -------------------------------------------------------------

# x<-quagum(runif(100000), c(1,3))
# x<-revd(10000,location= 1,scale= 3) # require(EnvStats)
# x<- rgumbel(10000, location=1,scale=3) # FAdist
# identical runs

#' @name Gumbel
#' @title Gumbel Distribution
#'
#' @description Two-parameter Gumbel (Extreme Value type I) distribution
#'   with location \eqn{\mu}{m} and scale \eqn{\sigma}{s}. The Gumbel
#'   distribution is the limiting distribution of block maxima and is
#'   widely used in extreme-value analysis of hydroclimatic variables.
#'   Density, distribution, quantile, and random generation functions are
#'   provided via the \pkg{FAdist} package. Fitting is performed via
#'   closed-form L-moments through \code{\link[lmom]{pelgum}}.
#'
#'   Define \eqn{z = (x - \mu)/\sigma}{z = (x-m)/s}. The probability
#'   density function is:
#'   \deqn{f(x) = \frac{1}{\sigma} \exp(-z - e^{-z})}{f(x) = exp(-z - exp(-z)) / s}
#'   where:
#'   * \eqn{\mu}{m} --- location parameter (\code{location})
#'   * \eqn{\sigma}{s} --- scale parameter (\code{scale})
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @export
#'

dgumbel=function(x,location,scale) {
  fx=FAdist::dgumbel(x=x,scale=scale,location=location)
  return(fx)
}

#' @rdname Gumbel
#' @export
pgumbel=function(q,location,scale) {
  FX=FAdist::pgumbel(q=q,scale=scale,location=location)
  return(FX)
}

#' @rdname Gumbel
#' @export
qgumbel=function(p,location,scale) {
  x=FAdist::qgumbel(p=p,scale=scale,location=location)
  return(x)
}

#' @rdname Gumbel
#' @export
rgumbel=function(n,location,scale) {
  x=FAdist::rgumbel(n=n,scale=scale,location=location)
  return(x)
}

#' @title fitlm_gumbel
#'
#' @description Fits the two-parameter Gumbel distribution to an xts
#'   series by the method of L-moments. Parameters are obtained in closed
#'   form via \code{\link[lmom]{pelgum}}, which matches the sample
#'   L-moment ratios to the theoretical L-moment ratios of the Gumbel
#'   distribution. Zero values below \code{zero_threshold} may be excluded
#'   via \code{ignore_zeros}. Goodness-of-fit is assessed via
#'   \code{GOF_tests}, and theoretical L-moments are computed from the
#'   fitted parameters for comparison with the sample.
#'
#' @param x An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted parameters), \code{TheorLMom} (theoretical L-moments),
#'   \code{DataLMom} (sample L-moments), and \code{GoF} (goodness-of-fit
#'   metrics).
#'
#' @examples
#' x <- xts::xts(rgumbel(365, location = 1, scale = 3),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' fit <- fitlm_gumbel(x)
#' fit$Param
#'
#' @export
#'

fitlm_gumbel=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- stats::na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }

  sam=lmom::samlmu(x)
  par=lmom::pelgum(sam)
  names(par) = NULL
  fit = list()
  fit$location=par[1]
  fit$scale = par[2]

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'gumbel')

  Res<-list()
  Res$Distribution<-list(FXs="qgumbel")
  Res$Param<-fit
  Res$TheorLMom<-lmom::lmrgum(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF<-GoF

  return(Res)
}



# LogNormal, 3-parameter -------------------------------------------------------------

# x<-qualn3(runif(10000), c(2,1,0.5))

# elmlnorm3(x,bound=0,plot=0)

#' @name Log-Normal
#' @title Log-Normal Distribution
#'
#' @description Three-parameter Log-Normal distribution with location
#'   \eqn{\mu}{m}, scale \eqn{\beta}{b}, and shape \eqn{\sigma}{s}. The
#'   distribution is lower-bounded at \eqn{\mu}{m} and arises when the
#'   logarithm of \eqn{x - \mu} follows a Normal distribution. Density,
#'   distribution, quantile, and random generation functions are provided
#'   via the \pkg{FAdist} package. Fitting is performed via closed-form
#'   L-moments through \code{\link[lmom]{pelln3}}; an optional fixed lower
#'   bound is supported.
#'
#'   The probability density function is:
#'   \deqn{f(x) = \frac{1}{(x-\mu) \sigma \sqrt{2\pi}} \exp\left(-\frac{\left(\ln\left(\frac{x-\mu}{\beta}\right)\right)^2}{2\sigma^2}\right), \quad x > \mu}{f(x) = exp(-(ln((x-m)/b))^2 / (2s^2)) / ((x-m) * s * sqrt(2*pi)), x > m}
#'   where:
#'   * \eqn{\mu}{m} --- location parameter (\code{location})
#'   * \eqn{\beta}{b} --- scale parameter (\code{scale})
#'   * \eqn{\sigma}{s} --- shape parameter (\code{shape})
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param shape shape parameter.
#' @export
#'

dlognorm=function(x,location=0,scale,shape) {
  fx=FAdist::dlnorm3(x=x,shape=shape,scale=scale,thres=location)
  return(fx)
}

#' @rdname Log-Normal
#' @export
plognorm=function(q,location=0,scale,shape) {
  FX=FAdist::plnorm3(q=q,shape=shape,scale=scale,thres=location)
  return(FX)
}

#' @rdname Log-Normal
#' @export
qlognorm=function(p,location=0,scale,shape) {
  x=FAdist::qlnorm3(p=p,shape=shape,scale=scale,thres=location)
  return(x)
}

#' @rdname Log-Normal
#' @export
rlognorm=function(n,location=0,scale,shape) {
  x=FAdist::rlnorm3(n=n,shape=shape,scale=scale,thres=location)
  return(x)
}

#' @title fitlm_lognorm
#'
#' @description Fits the three-parameter Log-Normal distribution to an
#'   xts series by the method of L-moments. Parameters are obtained in
#'   closed form via \code{\link[lmom]{pelln3}}, which matches the sample
#'   L-moment ratios to the theoretical L-moment ratios of the Log-Normal
#'   distribution. If \code{bound} is supplied the location is fixed;
#'   otherwise all three parameters are free. Zero values below
#'   \code{zero_threshold} may be excluded via \code{ignore_zeros}.
#'   Goodness-of-fit is assessed via \code{GOF_tests}.
#'
#' @param x An xts object containing the time series data.
#' @param bound Numeric or \code{NULL}. Optional fixed lower bound
#'   (location). Default \code{NULL}.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted parameters), \code{TheorLMom} (theoretical L-moments),
#'   \code{DataLMom} (sample L-moments), and \code{GoF} (goodness-of-fit
#'   metrics).
#'
#' @examples
#' x <- xts::xts(rlognorm(365, location = 0, scale = 1, shape = 0.5),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' fit <- fitlm_lognorm(x)
#' fit$Param
#'
#' @export
#'

fitlm_lognorm=function(x,bound=NULL,ignore_zeros = FALSE, zero_threshold = 0.01) {

  x <- stats::na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }

  sam=lmom::samlmu(x)
  par=lmom::pelln3(sam,bound=bound)
  names(par) = NULL
  fit = list()
  fit$location=par[1]
  fit$scale = par[2]
  fit$shape = par[3]

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'lognorm')

  Res<-list()
  Res$Distribution<-list(FXs="qlnorm",bound=bound)
  Res$Param<-fit
  Res$TheorLMom<-lmom::lmrln3(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF<-GoF

  return(Res)
}


# GEV -------------------------------------------------------------
# Need change - not ready

# x<-quagev(runif(100000), c(0,1,-0.5)) #lmom
# x<-rgev(...) # FAdist

#' @name GEV
#' @title Generalized Extreme Value Distribution
#'
#' @description Three-parameter Generalised Extreme Value (GEV)
#'   distribution with location \eqn{\mu}{m}, scale \eqn{\sigma}{s}, and
#'   shape \eqn{\xi}{j}. The GEV encompasses the Gumbel (\eqn{\xi = 0}),
#'   Frechet (\eqn{\xi > 0}), and reversed Weibull (\eqn{\xi < 0})
#'   families. Density, distribution, quantile, and random generation
#'   functions are implemented natively. Fitting is performed via
#'   closed-form L-moments through \code{\link[lmom]{pelgev}}.
#'
#'   Define \eqn{t(x) = (1 + \xi(x-\mu)/\sigma)^{-1/\xi}}{t(x) = (1 + j*(x-m)/s)^(-1/j)}. The probability density
#'   function is:
#'   \deqn{f(x) = \frac{1}{\sigma} t(x)^{\xi+1} \exp(-t(x)), \quad \xi \neq 0}{f(x) = t(x)^(j+1) * exp(-t(x)) / s}
#'   \deqn{f(x) = \frac{1}{\sigma} \exp\left(-\frac{x-\mu}{\sigma} - \exp\left(-\frac{x-\mu}{\sigma}\right)\right), \quad \xi = 0}{f(x) = exp(-(x-m)/s - exp(-(x-m)/s)) / s}
#'   with \eqn{1 + \xi(x-\mu)/\sigma > 0}{1 + j*(x-m)/s > 0}.
#'   where:
#'   * \eqn{\mu}{m} --- location parameter (\code{location})
#'   * \eqn{\sigma}{s} --- scale parameter (\code{scale})
#'   * \eqn{\xi}{j} --- shape parameter (\code{shape})
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param shape shape parameter.
#' @export
#'

dgev=function(x,location,scale,shape) {
  if (shape != 0){
    fx=1/scale*(1+shape*((x-location)/scale))^(-1/shape-1)*exp(-(1+shape*((x-location)/scale))^(-1/shape))
  }else{
    fx=exp(-(x - location)/scale)*exp(-exp(-(x - location)/scale))
  }

  return(fx)
}

#' @rdname GEV
#' @export
pgev=function(q,location,scale,shape) {
  if (shape != 0){
    FX= exp(-(1+shape*((q-location)/scale))^(-1/shape))
  }else{
    FX= exp(-exp(-(q - location)/scale))
  }
  return(FX)
}

#' @rdname GEV
#' @export
qgev = function(p,shape,scale,location){
    if (shape != 0){
      x <- location + scale/shape * ((-log(p))^(-shape) - 1)
    }else{
      x <- location - scale * log(-log(p))
    }

    return(x)
  }

#' @rdname GEV
#' @export
rgev=function(n,location,scale,shape) {
  x=qgev(stats::runif(n), location = location, scale = scale, shape = shape)
  return(x)
}

#' @title fitlm_gev
#'
#' @description Fits the three-parameter GEV distribution to an xts series
#'   by the method of L-moments. Parameters are obtained in closed form
#'   via \code{\link[lmom]{pelgev}}, which matches the sample L-moment
#'   ratios to the theoretical L-moment ratios of the GEV distribution.
#'   Zero values below \code{zero_threshold} may be excluded via
#'   \code{ignore_zeros}. Goodness-of-fit is assessed via
#'   \code{GOF_tests}, and theoretical L-moments are computed from the
#'   fitted parameters for comparison with the sample.
#'
#' @param x An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted parameters), \code{TheorLMom} (theoretical L-moments),
#'   \code{DataLMom} (sample L-moments), and \code{GoF} (goodness-of-fit
#'   metrics).
#'
#' @examples
#' x <- xts::xts(rgev(365, location = 0, scale = 1, shape = -0.5),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' fit <- fitlm_gev(x)
#' fit$Param
#'
#' @export
#'



fitlm_gev=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- stats::na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }

  sam=lmom::samlmu(x)
  par=lmom::pelgev(sam)
  names(par) = NULL
  fit = list()
  fit$location=par[1]
  fit$scale = par[2]
  fit$shape = par[3]

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'gev')

  Res<-list()
  Res$Distribution<-list(FXs="qgev")
  Res$Param<-fit
  Res$TheorLMom<-lmom::lmrgev(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF = GoF

  return(Res)
}



# GenPareto (need change) -------------------------------------------------------------


# lower bound at location parameter (false in lmom tutorial)

#' @name GenPareto
#' @title Generalized Pareto Distribution
#'
#' @description Three-parameter Generalised Pareto distribution (GPD) with
#'   location \eqn{\mu}{m}, scale \eqn{\sigma}{s}, and shape
#'   \eqn{\xi}{j}. The GPD is the limiting distribution of peaks over a
#'   threshold and is widely used in extreme-value modelling. The
#'   distribution is lower-bounded at the location parameter. Density,
#'   distribution, quantile, and random generation functions are provided
#'   via the \pkg{lmom} package. Fitting is performed via closed-form
#'   L-moments through \code{\link[lmom]{pelgpa}}; an optional fixed lower
#'   bound is supported.
#'
#'   Define \eqn{z = (x - \mu)/\sigma}{z = (x-m)/s}. The probability
#'   density function is:
#'   \deqn{f(x) = \frac{1}{\sigma} (1 - \xi z)^{1/\xi - 1}, \quad \xi \neq 0}{f(x) = (1 - j*z)^(1/j - 1) / s}
#'   \deqn{f(x) = \frac{1}{\sigma} e^{-z}, \quad \xi = 0}{f(x) = exp(-z) / s}
#'   with \eqn{z \ge 0}{z >= 0} and \eqn{1 - \xi z > 0}{1 - j*z > 0}.
#'   where:
#'   * \eqn{\mu}{m} --- location parameter (\code{location})
#'   * \eqn{\sigma}{s} --- scale parameter (\code{scale})
#'   * \eqn{\xi}{j} --- shape parameter (\code{shape})
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param shape shape parameter.
#' @export
#'

dgpd=function(x,location,scale,shape){
  if (shape != 0){
    fx = (1/scale)*(1 - shape*(x-location)/scale)^(1/shape - 1)
  } else {
    fx = (1/scale)*exp(-(x-location)/scale)
  }
  return(fx)
}

#' @rdname GenPareto
#' @export
pgpd=function(q,location,scale,shape){
  FX=lmom::cdfgpa(q, para=c(location,scale,shape))
  return(FX)
}

#' @rdname GenPareto
#' @export
qgpd=function(p,location,scale,shape){
  x=lmom::quagpa(p, para=c(location,scale,shape))
  return(x)
}

#' @rdname GenPareto
#' @export
rgpd=function(n,location,scale,shape){
  x=lmom::quagpa(stats::runif(n), para=c(location,scale,shape))
  return(x)
}

#' @title fitlm_GPD
#'
#' @description Fits the three-parameter Generalised Pareto distribution
#'   to an xts series by the method of L-moments. Parameters are obtained
#'   in closed form via \code{\link[lmom]{pelgpa}}, which matches the
#'   sample L-moment ratios to the theoretical L-moment ratios of the GPD.
#'   If \code{bound} is supplied the location is fixed; otherwise all three
#'   parameters are free. The shape parameter governs tail behaviour:
#'   values greater than -1/2 are required for finite variance, greater
#'   than -1/3 for finite skewness, and greater than -1/4 for finite
#'   kurtosis. Zero values below \code{zero_threshold} may be excluded via
#'   \code{ignore_zeros}. Goodness-of-fit is assessed via
#'   \code{GOF_tests}.
#'
#' @param x An xts object containing the time series data.
#' @param bound Numeric or \code{NULL}. Optional fixed lower bound
#'   (location). Default \code{NULL}.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted parameters), \code{TheorLMom} (theoretical L-moments),
#'   \code{DataLMom} (sample L-moments), and \code{GoF} (goodness-of-fit
#'   metrics).
#'
#' @examples
#' x <- xts::xts(rgpd(365, location = 0, scale = 1, shape = -0.2),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' fit <- fitlm_GPD(x)
#' fit$Param
#'
#' @export
#'

fitlm_GPD=function(x,bound=NULL,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- stats::na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }

  sam=lmom::samlmu(x)
  par=lmom::pelgpa(sam,bound=bound)
  names(par) = NULL
  fit = list()
  fit$location=par[1]
  fit$scale = par[2]
  fit$shape = par[3]

  print("Shape parameter must be >-1/2 for finite variance")
  print("Shape parameter must be >-1/3 for finite skewness")
  print("Shape parameter must be >-1/4 for finite kurtosis")
  
  GoF <- GOF_tests(x = x, fit = fit, distribution = 'gpd')

  Res<-list()
  Res$Distribution<-list(FXs="qgpd")
  Res$Param<-fit
  Res$TheorLMom<-lmom::lmrgpa(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF = GoF

  return(Res)
}


# Generalised Gamma -------------------------------------------------------------

#' @name GenGamma
#' @title Generalized Gamma Distribution
#'
#' @description Three-parameter Generalised Gamma (Stacy) distribution
#'   with scale \eqn{s}{s}, first shape \eqn{\alpha}{a} (\code{shape1}),
#'   and second shape \eqn{\beta}{b} (\code{shape2}). The distribution is
#'   supported on \eqn{x > 0} and contains the ordinary Gamma
#'   (\eqn{\beta = 1}) and Weibull (\eqn{\beta = \alpha}) as special
#'   cases. Density, distribution, quantile, and random generation
#'   functions are provided via \pkg{VGAM} using the Stacy
#'   parameterisation. Fitting is performed numerically via L-BFGS-B
#'   minimisation of a normalised L-moment error function, seeded from a
#'   pre-computed lookup table (\code{GG_InitValues}) with Gauss-Legendre
#'   quadrature evaluation of the L-moments at each optimisation step.
#'
#'   The probability density function is:
#'   \deqn{f(x) = \frac{\beta}{\Gamma(\alpha/\beta) \, s^\alpha} x^{\alpha-1} \exp\left(-\left(\frac{x}{s}\right)^\beta\right), \quad x > 0}{f(x) = b * x^(a-1) * exp(-(x/s)^b) / (Gamma(a/b) * s^a), x > 0}
#'   where:
#'   * \eqn{s}{s} --- scale parameter (\code{scale})
#'   * \eqn{\alpha}{a} --- first shape parameter (\code{shape1}); \eqn{k = \alpha/\beta}{k = a/b} is the Stacy shape
#'   * \eqn{\beta}{b} --- second shape parameter (\code{shape2}); the distribution reduces to Gamma(\eqn{\alpha}{a}, \eqn{s}{s}) when \eqn{\beta=1}{b=1}, and to Weibull(\eqn{s}{s}, \eqn{\alpha}{a}) when \eqn{\beta=\alpha}{b=a}
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param scale scale parameter.
#' @param shape1 first shape parameter.
#' @param shape2 second shape parameter.
#' @export
#'

dgengamma=function(x,scale, shape1, shape2){
  fx = VGAM::dgengamma.stacy(x = x,scale = scale, k = (shape1/shape2), d = shape2)
  return(fx)
}

#' @rdname GenGamma
#' @export
pgengamma=function(q,scale, shape1, shape2){
  FX = VGAM::pgengamma.stacy(q = q ,scale = scale, k = (shape1/shape2), d = shape2)
  return(FX)
}

#' @rdname GenGamma
#' @export
qgengamma=function(p,scale, shape1, shape2){
  X = VGAM::qgengamma.stacy(p =  p ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

#' @rdname GenGamma
#' @export
rgengamma=function(n,scale, shape1, shape2){
  X = VGAM::rgengamma.stacy(n =  n ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

#' @title fitlm_gengamma
#'
#' @description Fits the three-parameter Generalised Gamma (Stacy)
#'   distribution to an xts series by numerical minimisation of
#'   L-moments. As no closed-form L-moment expressions exist for this
#'   distribution, the L-moments are computed at each optimisation step
#'   via Gauss-Legendre quadrature of the PWMs. The L-BFGS-B optimiser
#'   minimises the normalised root-sum-square error between the sample and
#'   theoretical L-moments of orders specified by \code{order}. A
#'   two-step seeding procedure is employed: the shape parameters are
#'   initialised by a min-max nearest-neighbour search over the scale-free
#'   L-ratios in the \code{GG_InitValues} lookup table, and the scale is
#'   then derived analytically from the sample first L-moment and the
#'   matched shapes' unit-scale L-moment. Zero values below
#'   \code{zero_threshold} may be excluded via \code{ignore_zeros}.
#'
#' @param x An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#' @param order Integer vector of L-moment orders matched by the
#'   optimiser. Default \code{1:5}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted \code{scale}, \code{shape1}, \code{shape2}),
#'   \code{TheorLMom} (theoretical L-moments), \code{DataLMom} (sample
#'   L-moments), and \code{GoF} (goodness-of-fit metrics).
#'
#' @examples
#' x <- xts::xts(rgengamma(365, scale = 2, shape1 = 0.8, shape2 = 0.5),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' \dontrun{
#' fit <- fitlm_gengamma(x)
#' fit$Param
#' }
#'
#' @export
#'

fitlm_gengamma=function(x,ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:5))  {
  max_order = max(4,order)
  x <- stats::na.omit(coredata(x))
  PW = 1
  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }

  pfunction=pgengamma
  qfunction=qgengamma
  dfunction = dgengamma

  sample_LM = Lmoments::Lcoefs(NZ, rmax = max_order, na.rm = FALSE, trim = c(0, 0))

  init_values = GG_InitValues
  # Two-step seed (same as fitlm_dagum). (lcv, tau_3, tau_4) depend only on the shapes (scale
  # cancels in every ratio), so match the shapes by min-max NN on those scale-free quantities
  # over the unit-scale GG_InitValues, then derive scale = sample_lambda_1 / F, F the matched
  # shapes' unit-scale lambda_1. Avoids the blown-up absolute L-moments of the old table.
  lm_cols = c("lambda_1", "lambda_2", "tau_3", "tau_4", "lcv")
  rng = vapply(init_values[lm_cols], function(z) diff(range(z)), numeric(1))
  d = as.numeric(sample_LM)[1:4]
  d[5] <- d[2] / d[1]
  nn = which.min(((init_values$lcv - d[5]) / rng[5])^2 +
                 ((init_values$tau_3    - d[3]) / rng[3])^2 +
                 ((init_values$tau_4    - d[4]) / rng[4])^2)
  start_par = init_values[nn, c("scale", "shape1", "shape2")]
  unit_lms = lmom_gengamma(1:max(order), 1, as.numeric(start_par[2]), as.numeric(start_par[3]))
  start_par["scale"] = d[1] / unit_lms[1]
  start_par["scale"] = min(max(start_par["scale"], 0.001), 500)   # keep start in the optim box
  start_par = as.numeric(start_par)
  params_optim <- function(params, target_LMs, order){
    temp_lms <- as.numeric(lmom_gengamma(1:max(order), params[1], params[2], params[3]))
    # Degenerate shapes yield non-finite L-moments: large finite penalty so L-BFGS-B retreats.
    if (!all(is.finite(temp_lms))) return(1e6)
    temp_err <- sapply(order, FUN = function(x){((target_LMs[x] - temp_lms[x])/target_LMs[x])^2})
    temp_err <- sqrt(sum(temp_err))
    if (!is.finite(temp_err)) 1e6 else temp_err
  }
  all = stats::optim(start_par, fn = params_optim, target_LMs = sample_LM, order = order,
              method = "L-BFGS-B", lower = c(0.001, 0.005, 0.05), upper = c(500,2000,200))
  params = list(scale = all$par[1], shape1 = all$par[2], shape2 = all$par[3])

  TheorLmom = lmom_gengamma(1:5, scale=params$scale, shape1=params$shape1, shape2=params$shape2)

  GoF <- GOF_tests(x = NZ, fit = params, distribution = 'gengamma')

  #para$PW=PW

  Res<-list()
  Res$Distribution<-list(FXs="qgengamma")
  Res$Param<-params
  Res$TheorLMom<-TheorLmom
  Res$DataLMom<-as.data.frame(sample_LM)
  Res$GoF<-GoF
  #Res$Diag = c(itmax,NP)

  return(Res)
}


#Generalized Gamma with location----------------------------------------
#' @name GenGamma-Location
#' @title Generalized Gamma with Location Distribution
#'
#' @description Four-parameter Generalised Gamma distribution with
#'   location \eqn{\mu}{m}, scale \eqn{s}{s}, first shape
#'   \eqn{\alpha}{a}, and second shape \eqn{\beta}{b}. The distribution
#'   is lower-bounded at the location parameter and extends the
#'   three-parameter Stacy form by a shift. Fitting is delegated to
#'   \code{\link{fitlm_gengamma}} after subtracting the location from the
#'   data; the location is supplied by the user rather than fitted
#'   automatically. Density, distribution, quantile, and random generation
#'   functions are provided via \pkg{VGAM}.
#'
#'   The probability density function is:
#'   \deqn{f(x) = \frac{\beta}{\Gamma(\alpha/\beta) \, s^\alpha} (x-\mu)^{\alpha-1} \exp\left(-\left(\frac{x-\mu}{s}\right)^\beta\right), \quad x \ge \mu}{f(x) = b * (x-m)^(a-1) * exp(-((x-m)/s)^b) / (Gamma(a/b) * s^a), x >= m}
#'   where:
#'   * \eqn{\mu}{m} --- location parameter (\code{location})
#'   * \eqn{s}{s} --- scale parameter (\code{scale})
#'   * \eqn{\alpha}{a} --- first shape parameter (\code{shape1})
#'   * \eqn{\beta}{b} --- second shape parameter (\code{shape2})
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param shape1 first shape parameter.
#' @param shape2 second shape parameter.
#' @export
#'

pgengamma_loc=function(q, location, scale, shape1, shape2){
  FX = VGAM::pgengamma.stacy(q = q-location, scale = scale, k = (shape1/shape2), d = shape2)
  return(FX)
}

#' @rdname GenGamma-Location
#' @export
dgengamma_loc=function(x,location, scale, shape1, shape2){
  fx = VGAM::dgengamma.stacy(x = x - location,scale = scale, k = (shape1/shape2), d = shape2)
  return(fx)
}

#' @rdname GenGamma-Location
#' @export
qgengamma_loc=function(p, location, scale, shape1, shape2){
  X = location+ VGAM::qgengamma.stacy(p =  p ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

#' @rdname GenGamma-Location
#' @export
rgengamma_loc=function(n, location, scale, shape1, shape2){
  X = location + VGAM::rgengamma.stacy(n =  n ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

#' @title fitlm_gengamma_loc
#'
#' @description Fits the four-parameter Generalised Gamma distribution
#'   with location to an xts series. The location parameter is supplied by
#'   the user and is subtracted from the data before fitting; the
#'   remaining three parameters (scale, shape1, shape2) are then fitted by
#'   delegating to \code{\link{fitlm_gengamma}}. The theoretical first
#'   L-moment is adjusted by adding back the location after the fit. Zero
#'   values below \code{zero_threshold} may be excluded via
#'   \code{ignore_zeros}. Goodness-of-fit is assessed via
#'   \code{GOF_tests} evaluated on the full four-parameter set.
#'
#' @param x An xts object containing the time series data.
#' @param location Numeric. The location parameter of the distribution.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#' @param order Integer vector of L-moment orders passed to
#'   \code{fitlm_gengamma}. Default \code{1:5}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted \code{location}, \code{scale}, \code{shape1},
#'   \code{shape2}), \code{TheorLMom} (theoretical L-moments),
#'   \code{DataLMom} (sample L-moments), and \code{GoF} (goodness-of-fit
#'   metrics).
#'
#' @examples
#' x <- xts::xts(rgengamma_loc(365, location = 1, scale = 2,
#'         shape1 = 0.8, shape2 = 0.5),
#'         order.by = as.Date("2020-01-01") + 0:364)
#' \dontrun{
#' fit <- fitlm_gengamma_loc(x, location = 1)
#' fit$Param
#' }
#'
#' @export
#'


fitlm_gengamma_loc=function(x, location, ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:5))  {
  max_order = max(4,order)
  x <- stats::na.omit(coredata(x))
  PW = 1
  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }

  sample_LM = Lmoments::Lcoefs(NZ, rmax = max_order, na.rm = FALSE, trim = c(0, 0))

  # pass a 1-column matrix so the recursive call's `x[x > zt, ]` indexing stays 2-D
  temp_fit = fitlm_gengamma(x = matrix(NZ - location, ncol = 1), ignore_zeros = ignore_zeros, zero_threshold = zero_threshold, order = order)

  TheorLmom= temp_fit$TheorLMom
  TheorLmom[1] = TheorLmom[1] + location

  full_par = c(list(location = location), temp_fit$Param)

  GoF <- GOF_tests(x = NZ, fit = full_par, distribution = 'gengamma_loc')

  #para$PW=PW

  Res<-list()
  Res$Distribution<-list(FXs="qgengamma_loc")
  Res$Param<-full_par
  Res$TheorLMom<-TheorLmom
  Res$DataLMom<-as.data.frame(sample_LM)
  Res$GoF<-GoF
  #Res$Diag = c(itmax,NP)

  return(Res)
}



# Burr-XII -------------------------------------------------------------

# dburr=function(x,scale,shape1,shape2) {
#   fx=ExtDist::dBurr(x = x,b=scale,g=shape1,s=shape2)
#   return(fx)
# }
#
# pburr=function(q,scale,shape1,shape2) {
#   FX=ExtDist::pBurr(q = q,b=scale,g=shape1,s=shape2)
#   return(FX)
# }
#
# qburr=function(p,scale,shape1,shape2) {
#   x=ExtDist::qBurr(p = p,b=scale,g=shape1,s=shape2)
#   return(x)
#   # u=0.1;
#   # T=1/(1-(u));
#   # scale=1.5; shape1=2; shape2=3;
#   # qburr(u, scale, shape1, shape2);
#   # scale*(T^(1/shape2) -1)^(1/shape1);
#   # scale*((1-u)^(-1/shape2) -1)^(1/shape1)
# }
# rburr=function(n,scale,shape1,shape2) {
#   x=ExtDist::rBurr(n=n,b=scale,g=shape1,s=shape2)
#   return(x)
# }

#DK definition

#' @name BurXII
#' @title Burr Type XII Distribution
#'
#' @description Three-parameter Burr Type XII distribution with scale
#'   \eqn{s}{s}, first shape \eqn{\zeta}{z} (\code{shape1}), and
#'   second shape \eqn{\theta}{th} (\code{shape2}, the tail index). The
#'   distribution is supported on \eqn{x > 0} and is widely used for
#'   modelling heavy-tailed positive variables. L-moments are computed via
#'   closed-form expressions involving Beta-function PWMs; a two-step
#'   seeding procedure initialises the shapes from a pre-computed lookup
#'   table (\code{Burr_InitValues}) and derives the scale analytically.
#'   Fitting uses L-BFGS-B minimisation of the normalised L-moment error.
#'
#'   The probability density function is:
#'   \deqn{f(x) = \zeta s^{-\zeta} x^{\zeta-1} \left(\zeta\theta\left(\frac{x}{s}\right)^\zeta + 1\right)^{-1/(\zeta\theta) - 1}, \quad x > 0}{f(x) = z*s^(-z) * x^(z-1) * (z*th*(x/s)^z + 1)^(-1/(z*th) - 1), x > 0}
#'   where:
#'   * \eqn{s}{s} --- scale parameter (\code{scale})
#'   * \eqn{\zeta}{z} --- first shape parameter (\code{shape1}), controls lower-tail behaviour
#'   * \eqn{\theta}{th} --- second shape parameter (\code{shape2}), tail index; the mean exists only when \eqn{\theta < 1}{th < 1}
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param scale scale parameter.
#' @param shape1 first shape parameter.
#' @param shape2 second shape parameter.
#' @param PW probability weight (point mass at zero for zero-inflated variants).
#' @export
#'

dburr=function(x, scale, shape1, shape2, PW=1){
  d=PW*shape1*scale^(-shape1)*x^(shape1-1)*((shape1*shape2*(x/scale)^shape1)+1)^(-1/(shape1*shape2)-1)
}

#' @rdname BurXII
#' @export
pburr=function(q, scale, shape1, shape2, PW=1) {
  # scale: lamda
  # shape1: zeta
  # shape2: tail index!
  p=1-PW*(1+(shape1*shape2)*(q/scale)^shape1)^(-1/(shape1*shape2))
  return(p)
}

#' @rdname BurXII
#' @export
qburr=function(p, scale, shape1, shape2, PW=1) {

  q=rep(NA, length(p))

  for (i in 1:length(p) ) {
    if (p[i]>(1-PW)) {
      A=1-((1-p[i])/PW)^(-(shape1*shape2))
      q[i]=scale*(-A/(shape1*shape2))^(1/shape1)
    } else {
      q[i]=0
    }
  }
  return(q)
}

#' @rdname BurXII
#' @export
rburr=function(n, scale, shape1, shape2, PW=1) {
  p=stats::runif(n)
  q=qburr(p = p, scale=scale, shape1=shape1, shape2=shape2, PW=PW)
  return(q)
}

#' @title fitlm_burr
#'
#' @description Fits the three-parameter Burr Type XII distribution to an
#'   xts series by numerical minimisation of L-moments. Closed-form
#'   L-moment expressions via Beta-function PWMs are evaluated at each
#'   optimisation step. The L-BFGS-B optimiser minimises the normalised
#'   root-sum-square error between the sample and theoretical L-moments of
#'   orders specified by \code{order}. A two-step seeding procedure
#'   initialises the shape parameters by min-max nearest-neighbour search
#'   over the scale-free L-ratios in the \code{Burr_InitValues} lookup
#'   table, and the scale is then derived analytically from the sample
#'   first L-moment and the matched shapes' unit-scale L-moment. The
#'   mean exists only when \eqn{\theta < 1}, imposing an effective upper
#'   bound on the tail index during optimisation. Zero values below
#'   \code{zero_threshold} may be excluded via \code{ignore_zeros}.
#'
#' @param x An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#' @param order Integer vector of L-moment orders matched by the
#'   optimiser. Default \code{1:5}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted \code{scale}, \code{shape1}, \code{shape2}),
#'   \code{TheorLMom} (theoretical L-moments), \code{DataLMom} (sample
#'   L-moments), and \code{GoF} (goodness-of-fit metrics).
#'
#' @examples
#' x <- xts::xts(rburr(365, scale = 1.5, shape1 = 2, shape2 = 0.3),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' \dontrun{
#' fit <- fitlm_burr(x)
#' fit$Param
#' }
#'
#' @export
#'

fitlm_burr=function(x,ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:5))  {
  max_order = max(4,order)
  x <- stats::na.omit(zoo::coredata(x))
  PW = 1
  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }

  sample_LM = Lmoments::Lcoefs(NZ, rmax = max_order, na.rm = FALSE, trim = c(0, 0))
  
  init_values = Burr_InitValues
  # Two-step seed (same as fitlm_dagum). (lcv, tau_3, tau_4) depend only on the shapes (scale
  # cancels in every ratio), so match the shapes by min-max NN on those scale-free quantities
  # over the unit-scale Burr_InitValues, then derive scale = sample_lambda_1 / F, F the matched
  # shapes' unit-scale lambda_1 -- avoiding the earlier four-column metric that mixed the
  # scale-carrying lambda_1 with the scale-free ratios.
  lm_cols = c("lambda_1", "lambda_2", "tau_3", "tau_4", "lcv")
  rng = vapply(init_values[lm_cols], function(z) diff(range(z)), numeric(1))
  d = as.numeric(sample_LM)[1:4]
  d[5] <- d[2] / d[1]
  nn = which.min(((init_values$lcv - d[5]) / rng[5])^2 +
                 ((init_values$tau_3    - d[3]) / rng[3])^2 +
                 ((init_values$tau_4    - d[4]) / rng[4])^2)
  start_par = init_values[nn, c("scale", "shape1", "shape2")]
  unit_lms = lmom_burr(1:max(order), 1, as.numeric(start_par[2]), as.numeric(start_par[3]))
  start_par["scale"] = d[1] / unit_lms[1]
  start_par["scale"] = min(max(start_par["scale"], 0.5), 50)   # keep start in the optim box
  start_par = as.numeric(start_par)
  params_optim <- function(params, target_LMs, order){
    temp_lms <- as.numeric(lmom_burr(1:max(order), params[1], params[2], params[3]))
    # Out-of-domain (shape2 -> 1, mean diverges): closed form yields Inf. Return a large
    # finite penalty so L-BFGS-B retreats from the edge instead of failing on a non-finite fn.
    if (!all(is.finite(temp_lms))) return(1e6)
    temp_err <- sapply(order, FUN = function(x){((target_LMs[x] - temp_lms[x])/target_LMs[x])^2})
    temp_err <- sqrt(sum(temp_err))
    if (!is.finite(temp_err)) 1e6 else temp_err
  }
  all = stats::optim(start_par, fn = params_optim, target_LMs = sample_LM, order = order,
              method = "L-BFGS-B", lower = c(0.5, 0.5, 0.001), upper = c(50,50,1))
  params = list(scale = all$par[1], shape1 = all$par[2], shape2 = all$par[3])

  TheorLmom=lmom_burr(1:5, scale=params$scale, shape1=params$shape1, shape2=params$shape2)

  GoF <- GOF_tests(x = NZ, fit = params, distribution = 'burr')

  #para$PW=PW

  Res<-list()
  Res$Distribution<-list(FXs="qburr")
  Res$Param<-params
  Res$TheorLMom<-TheorLmom
  Res$DataLMom<-sample_LM
  Res$GoF<-GoF
  #Res$Diag = c(itmax,NP)

  return(Res)
}


#Dagum--------------------------------------------------------------

#' @name Dagum
#' @title Dagum Distribution
#'
#' @description Three-parameter Dagum (Burr Type III) distribution with
#'   scale \eqn{s}{s}, first shape \eqn{\alpha}{a} (\code{shape1}), and
#'   second shape \eqn{\theta}{th} (\code{shape2}, the tail index). The
#'   distribution is supported on \eqn{x > 0} and, unlike the Burr XII,
#'   has a unimodal density with mode at the origin for certain parameter
#'   ranges. L-moments are computed via closed-form Beta-function PWMs;
#'   fitting uses L-BFGS-B minimisation seeded from the
#'   \code{Dagum_InitValues} lookup table with scale-free L-ratio matching
#'   and analytic scale derivation.
#'
#'   The probability density function is:
#'   \deqn{f(x) = \frac{1}{s} \left(\frac{x}{s}\right)^{-1/\theta - 1} \left(1 + \frac{\theta}{\alpha}\left(\frac{x}{s}\right)^{-1/\theta}\right)^{-\alpha-1}, \quad x > 0}{f(x) = 1/s * (x/s)^(-1/th - 1) * (1 + th/a * (x/s)^(-1/th))^(-a-1), x > 0}
#'   where:
#'   * \eqn{s}{s} --- scale parameter (\code{scale})
#'   * \eqn{\alpha}{a} --- first shape parameter (\code{shape1}), controls upper-tail heaviness
#'   * \eqn{\theta}{th} --- second shape parameter (\code{shape2}), tail index; the mean exists only when \eqn{\theta < 1}{th < 1}
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param scale scale parameter.
#' @param shape1 first shape parameter.
#' @param shape2 second shape parameter.
#' @param PW probability weight (point mass at zero for zero-inflated variants).
#' @export
#'

pdagum=function(q, scale, shape1, shape2, PW=1) {
  # The moments exist (i.e., have finite values) only for order r < 1/shape2;  for larger r they are infinite.
  p=(1-PW)+PW*(1+(shape2/shape1)*(q/scale)^(-1/shape2))^-shape1
  return(p)
}

#' @rdname Dagum
#' @export
ddagum=function(x, scale, shape1, shape2, PW=1) {
  d=PW*(1/scale)*(x/scale)^(-1-1/shape2)*((1/shape1)*(shape1 + shape2*(x/scale)^(-1/shape2)))^(-1-shape1)
}

#' @rdname Dagum
#' @export
qdagum=function(p, scale, shape1, shape2, PW=1) {

  q=rep(NA, length(p))
  q=scale * (-(((1 - (-((1 - PW - p)/PW))^(-1/shape1)) *shape1)/shape2))^-shape2

  for (i in 1:length(p) ) {
    if (p[i]>(1-PW)) {
      q[i]=scale * (-(((1 - (-((1 - PW - p[i])/PW))^(-1/shape1)) *shape1)/shape2))^-shape2
    } else {
      q[i]=0
    }
  }
  return(q)
}

#' @rdname Dagum
#' @export
rdagum=function(n, scale, shape1, shape2, PW=1) {
  qdagum(p = stats::runif(n), scale=scale, shape1 = shape1, shape2 = shape2, PW=PW)
}


#' @title fitlm_dagum
#'
#' @description Fits the three-parameter Dagum distribution to an xts
#'   series by numerical minimisation of L-moments. Closed-form L-moment
#'   expressions via Beta-function PWMs are evaluated at each optimisation
#'   step. The L-BFGS-B optimiser minimises the normalised root-sum-square
#'   error between the sample and theoretical L-moments of orders
#'   specified by \code{order}. A two-step seeding procedure initialises
#'   the shape parameters by min-max nearest-neighbour search over the
#'   scale-free L-ratios in the \code{Dagum_InitValues} lookup table, and
#'   the scale is derived analytically from the sample first L-moment and
#'   the matched shapes' unit-scale L-moment. The mean exists only when
#'   \eqn{\theta < 1}. Zero values below \code{zero_threshold} may be
#'   excluded via \code{ignore_zeros}.
#'
#' @param x An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#' @param order Integer vector of L-moment orders matched by the
#'   optimiser. Default \code{1:5}.
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted \code{scale}, \code{shape1}, \code{shape2}),
#'   \code{TheorLMom} (theoretical L-moments), \code{DataLMom} (sample
#'   L-moments), and \code{GoF} (goodness-of-fit metrics).
#'
#' @examples
#' x <- xts::xts(rdagum(365, scale = 1.5, shape1 = 2, shape2 = 0.3),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' \dontrun{
#' fit <- fitlm_dagum(x)
#' fit$Param
#' }
#'
#' @export
#'

fitlm_dagum=function(x,ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:5))  {
  max_order = max(4,order)
  x <- stats::na.omit(coredata(x))

  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }

  sample_LM = Lmoments::Lcoefs(NZ, rmax = max_order, na.rm = FALSE, trim = c(0, 0))
  
  init_values = Dagum_InitValues
  # Seed optim from the init row nearest in the four absolute L-moments
  # (lambda_1, lambda_2, tau_3, tau_4), min-max scaled so each column contributes comparably.
  # Min-max scaling reduces to weighting each squared term by 1/range^2 (the minima cancel
  # in the squared distance), so only the column ranges are needed. Because lambda_1 (the
  # scale carrier) is part of the metric, the matched row's parameters seed optim directly,
  # by name, with no scale derivation -- avoiding the earlier scale-free-ratio NN that left
  # optim with an arbitrary scale and trapped shape1 >> 100.
  lm_cols = c("lambda_1", "lambda_2", "tau_3", "tau_4", "lcv")
  rng = vapply(init_values[lm_cols], function(z) diff(range(z)), numeric(1))
  d = as.numeric(sample_LM)[1:4]
  d[5] <- d[2] / d[1]
  nn = which.min(((init_values$lcv - d[5]) / rng[5])^2 +
                 ((init_values$tau_3    - d[3]) / rng[3])^2 +
                 ((init_values$tau_4    - d[4]) / rng[4])^2)
  # nn = which.min(((init_values$lambda_1 - d[1]) / rng[1])^2 +
  #                ((init_values$lambda_2 - d[2]) / rng[2])^2 +
  #                ((init_values$tau_3    - d[3]) / rng[3])^2 +
  #                ((init_values$tau_4    - d[4]) / rng[4])^2)
  start_par = init_values[nn, c("scale", "shape1", "shape2")]
  unit_lms = lmom_dagum(1:max(order), 1, as.numeric(start_par[2]), as.numeric(start_par[3]))
  start_par["scale"] = d[1] / unit_lms[1]
  start_par["scale"] = min(max(start_par["scale"], 0.5), 1000)   # keep start in the optim box
  start_par = as.numeric(start_par)
  params_optim <- function(params, target_LMs, order){
    temp_lms <- as.numeric(lmom_dagum(1:max(order), params[1], params[2], params[3]))
    # Out-of-domain (shape2 -> 1, mean diverges): closed form yields Inf. Return a large
    # finite penalty so L-BFGS-B retreats from the edge instead of failing on a non-finite fn.
    if (!all(is.finite(temp_lms))) return(1e6)
    temp_err <- sapply(order, FUN = function(x){((target_LMs[x] - temp_lms[x])/target_LMs[x])^2})
    temp_err <- sqrt(sum(temp_err))
    if (!is.finite(temp_err)) 1e6 else temp_err
  }
  all = stats::optim(start_par, fn = params_optim, target_LMs = sample_LM, order = order,
              method = "L-BFGS-B", lower = c(0.5, 0.0001, 0.000001), upper = c(1000,500,0.5))
  params = list(scale = all$par[1], shape1 = all$par[2], shape2 = all$par[3])

  TheorLmom=lmom_dagum(1:5, scale=params$scale, shape1=params$shape1, shape2=params$shape2)

  GoF <- GOF_tests(x = NZ, fit = params, distribution = 'dagum')

  #para$PW=PW

  Res<-list()
  Res$Distribution<-list(FXs="qdagum")
  Res$Param<-params
  Res$TheorLMom<-TheorLmom
  Res$DataLMom<-as.data.frame(sample_LM)
  Res$GoF<-GoF
  #Res$Diag = c(itmax,NP)

  return(Res)
}


# Exponentiated Weibull (Generalized Weibull) Distribution-----------------------------
# Cumulative distribution function
#' @name ExpWeibull
#' @title Exponential Weibull Distribution
#'
#' @description Three-parameter Exponentiated Weibull distribution with
#'   scale \eqn{s}{s}, first shape \eqn{a}{a} (\code{shape1}, the
#'   Weibull shape), and second shape \eqn{b}{b} (\code{shape2}, the
#'   exponentiation power). The distribution is supported on \eqn{x > 0}
#'   and contains the ordinary Weibull (\eqn{b = 1}) and Exponential
#'   (\eqn{a = 1, b = 1}) as special cases. L-moments are computed by
#'   tanh-sinh (double-exponential) quadrature, which is well-conditioned
#'   across the full positive parameter range. Fitting uses L-BFGS-B
#'   minimisation seeded from the \code{ExpWeibull_InitValues} lookup
#'   table with scale-free L-ratio matching and analytic scale derivation.
#'
#'   The probability density function is:
#'   \deqn{f(x) = \frac{b a}{s} \left(\frac{x}{s}\right)^{a-1} \exp\left(-\left(\frac{x}{s}\right)^a\right) \left(1 - \exp\left(-\left(\frac{x}{s}\right)^a\right)\right)^{b-1}, \quad x > 0}{f(x) = b*a/s * (x/s)^(a-1) * exp(-(x/s)^a) * (1-exp(-(x/s)^a))^(b-1), x > 0}
#'   where:
#'   * \eqn{s}{s} --- scale parameter (\code{scale})
#'   * \eqn{a}{a} --- first shape parameter (\code{shape1}), the Weibull shape; a = 1 gives the exponential baseline
#'   * \eqn{b}{b} --- second shape parameter (\code{shape2}), the exponentiation power; b = 1 recovers the ordinary Weibull
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken as the number required.
#' @param scale scale parameter.
#' @param shape1 first shape parameter.
#' @param shape2 second shape parameter.
#' @param log,log.p logical; if \code{TRUE}, probabilities p are given as log(p).
#' @export
#'

pexpweibull<- function(q,scale,shape1,shape2,log.p=FALSE){
  log.cdf <- shape2*stats::pweibull(q = q, scale=scale,shape=shape1,log.p=TRUE)
  ifelse(log.p, return(log.cdf), return(exp(log.cdf)))
}

# Probability density function
#' @rdname ExpWeibull
#' @export
dexpweibull<- function(x,scale,shape1,shape2,log=FALSE){
  log.pdf <-  log(shape2) + (shape2-1)*stats::pweibull(q = x,scale=scale,shape=shape1,log=TRUE) +
    stats::dweibull(x = x,scale=scale,shape=shape1,log=TRUE)
  ifelse(log, return(log.pdf), return(exp(log.pdf)))
}

# Quantile function
#' @rdname ExpWeibull
#' @export
qexpweibull<- function(p,scale,shape1,shape2){
  quant <-  stats::qweibull(p = p^(1/shape2),scale=scale,shape=shape1)
  return(quant)
}

# Random number generation
#' @rdname ExpWeibull
#' @export
rexpweibull<- function(n,scale,shape1,shape2){
  u = stats::runif(n)
  sim <-  stats::qweibull(u^(1/shape2),scale=scale,shape=shape1)
  return(sim)
}

#' @title fitlm_expweibull
#'
#' @description Fits the three-parameter Exponentiated Weibull
#'   distribution to an xts series by numerical minimisation of
#'   L-moments. L-moments are computed at each optimisation step via a
#'   fast tanh-sinh (double-exponential) quadrature engine
#'   (\code{lmom_expweibull}), chosen for its well-conditioned behaviour
#'   across the full positive parameter range. The L-BFGS-B optimiser
#'   minimises the normalised root-sum-square error between the sample and
#'   theoretical L-moments of orders specified by \code{order}. A
#'   two-step seeding procedure initialises the shape parameters from the
#'   \code{ExpWeibull_InitValues} lookup table and derives the scale
#'   analytically. The distribution is meant to be fitted only to
#'   non-negative data. Zero values below \code{zero_threshold} may be
#'   excluded via \code{ignore_zeros}.
#'
#' @param x An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#' @param order Integer vector of L-moment orders matched by the
#'   optimiser. Default \code{1:3} (exact method-of-L-moments for the
#'   three parameters).
#'
#' @return A list with elements \code{Distribution}, \code{Param} (named
#'   list of fitted \code{scale}, \code{shape1}, \code{shape2}),
#'   \code{TheorLMom} (theoretical L-moments), \code{DataLMom} (sample
#'   L-moments), and \code{GoF} (goodness-of-fit metrics).
#'
#' @examples
#' x <- xts::xts(rexpweibull(365, scale = 3, shape1 = 1.5, shape2 = 2),
#'               order.by = as.Date("2020-01-01") + 0:364)
#' \dontrun{
#' fit <- fitlm_expweibull(x)
#' fit$Param
#' }
#'
#' @export
#'

fitlm_expweibull=function(x,ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:3)) {
  # Important: MEANT to be fitted ONLY to NON-NEGATIVE data.
  max_order = max(4,order)
  x <- stats::na.omit(coredata(x))

  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }

  sample_LM = Lmoments::Lcoefs(NZ, rmax = max_order, na.rm = FALSE, trim = c(0, 0))

  init_values = ExpWeibull_InitValues
  # Two-step seed (same as fitlm_gengamma / fitlm_dagum). (lcv, tau_3, tau_4) depend ONLY on the
  # shapes (a, b) -- scale cancels in every ratio -- so match the shapes by min-max NN on those
  # three scale-free quantities, then derive the scale: scale = sample_lambda_1 / F, where F is the
  # matched shapes' unit-scale lambda_1.
  lm_cols = c("lambda_1", "lambda_2", "tau_3", "tau_4", "lcv")
  rng = vapply(init_values[lm_cols], function(z) diff(range(z)), numeric(1))
  d = as.numeric(sample_LM)[1:4]
  d[5] <- d[2] / d[1]
  nn = which.min(((init_values$lcv - d[5]) / rng[5])^2 +
                 ((init_values$tau_3    - d[3]) / rng[3])^2 +
                 ((init_values$tau_4    - d[4]) / rng[4])^2)
  start_par = init_values[nn, c("scale", "shape1", "shape2")]
  unit_lms = lmom_expweibull(1:max(order), 1, as.numeric(start_par[2]), as.numeric(start_par[3]))
  start_par["scale"] = d[1] / unit_lms[1]
  start_par["scale"] = min(max(start_par["scale"], 0.001), 1e6)   # keep start in the optim box
  start_par = as.numeric(start_par)

  params_optim <- function(params, target_LMs, order){
    temp_lms <- as.numeric(lmom_expweibull(1:max(order), params[1], params[2], params[3]))
    if (!all(is.finite(temp_lms))) return(1e6)
    temp_err <- sapply(order, FUN = function(x){((target_LMs[x] - temp_lms[x])/target_LMs[x])^2})
    temp_err <- sqrt(sum(temp_err))
    if (!is.finite(temp_err)) 1e6 else temp_err
  }
  all = stats::optim(start_par, fn = params_optim, target_LMs = sample_LM, order = order,
              method = "L-BFGS-B", lower = c(0.001, 0.05, 0.02), upper = c(1e6,50,50))
  params = list(scale = all$par[1], shape1 = all$par[2], shape2 = all$par[3])

  TheorLmom=lmom_expweibull(1:5, scale=params$scale, shape1=params$shape1, shape2=params$shape2)

  GoF <- GOF_tests(x = NZ, fit = params, distribution = 'expweibull')

  Res<-list()
  Res$Distribution<-list(FXs="qexpweibull")
  Res$Param<-params
  Res$TheorLMom<-TheorLmom
  Res$DataLMom<-as.data.frame(sample_LM)
  Res$GoF<-GoF

  return(Res)

}
