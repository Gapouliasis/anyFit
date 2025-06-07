GOF_tests = function(x, fit, distribution){
  u_emp <- rank(coredata(x), na.last = NA, ties.method = "average")/(length(x)+1)
  q_emp <- coredata(x)
  qq_fitted <- do.call(paste0('q',distribution), c(list(p = u_emp), fit))

  CM <- CDFt::CramerVonMisesTwoSamples(q_emp, qq_fitted)
  KS <- CDFt::KolmogorovSmirnov(q_emp, qq_fitted)

  MLE = -sum(log(do.call(paste0('d',distribution),c(list(x = x), fit)) + 0.000001))
  plotpos<-lmomco::pp(x=x,a=0,sort=FALSE)
  theorquantiles = do.call(paste0('q',distribution), c(list(p = plotpos), fit))
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]

  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10

  # AIC = 2*length(fit)-2*sum(log(do.call(paste0('d',distribution),c(list(x = x), fit))))
  # BIC = length(fit)*log(length(x))-2*sum(do.call(paste0('d',distribution),c(list(x = x), fit)))

  GoF <- list(MLE=MLE, CM = CM, KS = KS,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  return(GoF)
}


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
#' @description Exponential distribution
#'
#' @import lmom
#'
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
  x=lmom::quaexp(f=runif(n),para=c(location,scale))
  return(x)
}


#' @title fitlm_exp
#'
#' @description Function for fitting the Exponential distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param bound Is the distribution bound? Default is NULL.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @export
#'

fitlm_exp=function(x,bound=NULL, ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- na.omit(coredata(x))
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
#' @description Rayleigh distribution
#'
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
#' @description Function for fitting the Rayleigh distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @export
#'

fitlm_rayleigh = function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- na.omit(coredata(x))
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
#' @description Gamma distribution
#'
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
#' @description Function for fitting the Gamma distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @export
#'

fitlm_gamma=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- na.omit(coredata(x))
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
  Res$TheorLMom<-lmrgam(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF<-GoF


  return(Res)
}


# Gamma3 -------------------------------------------------------------

# x<-rgamma3(10000,location=5,scale=0.5,shape=3)
#' @name Gamma3
#' @title 3-parameter Gamma Distribution
#'
#' @description 3-parameter Gamma distribution
#'
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
#' @description Function for fitting the 3-parameter Gamma distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @export
#'


fitlm_gamma3=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- na.omit(coredata(x))
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

#' @title fitlm_genlogi
#'
#' @description Function for fitting the Generalized Logistic distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @export
#'


fitlm_genlogi=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- na.omit(coredata(x))
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

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'gamma3')

  Res<-list()
  Res$Param<-fit
  Res$TheorLMom<-lmrglo(par,nmom=4)
  Res$DataLMom<-sam

  return(Res)
}


# Normal -------------------------------------------------------------

# x<-quanor(runif(10000), c(0,3))
# x<-rnorm(10000, mean=0,sd=3)
# identical runs

#' @name Normal
#' @title Normal Distribution
#'
#' @description Normal distribution
#'
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
#' @description Function for fitting the Normal distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @export
#'


fitlm_norm=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- na.omit(coredata(x))
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
#' @description 3-parameter Weibull distribution
#'
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
#' @description Function for fitting the 3-parameter Weibull distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param bound Is the distribution bound? Default is NULL.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @export
#'

fitlm_weibull=function(x,bound=NULL,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- na.omit(coredata(x))
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
  Res$TheorLMom<-lmrwei(par,nmom=4)
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
#' @description Gumbel distribution
#'
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
#' @description Function for fitting the Gumbel distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @export
#'

fitlm_gumbel=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- na.omit(coredata(x))
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
  Res$TheorLMom<-lmrgum(par,nmom=4)
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
#' @description Log-Normal distribution
#'
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
#' @description Function for fitting the 3-parameter Log-Normal distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param bound Is the distribution bound? Default is NULL.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @export
#'

fitlm_lognorm=function(x,bound=NULL,ignore_zeros = FALSE, zero_threshold = 0.01) {

  x <- na.omit(coredata(x))
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
  Res$TheorLMom<-lmrln3(par,nmom=4)
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
#' @description Generalized Extreme Value distribution
#'
#' @export
#'

dgev=function(x,location,scale,shape) {
  fx=FAdist::dgev(x=x,shape=shape,scale=scale,location=location)
  return(fx)
}

#' @rdname GEV
#' @export
pgev=function(q,location,scale,shape) {
  FX=FAdist::pgev(q=q,shape=shape,scale=scale,location=location)
  return(FX)
}

#' @rdname GEV
#' @export
qgev=function(p,location,scale,shape) {
  x=FAdist::qgev(p=p,shape=shape,scale=scale,location=location)
  return(x)
}

#' @rdname GEV
#' @export
rgev=function(n,location,scale,shape) {
  x=FAdist::rgev(n=n,shape=shape,scale=scale,location=location)
  return(x)
}

#' @title fitlm_gev
#'
#' @description Function for fitting the Generalized Extreme Value distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @export
#'


fitlm_gev=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- na.omit(coredata(x))
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
  Res$TheorLMom<-lmrgev(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF = GoF

  return(Res)
}



# GenPareto (need change) -------------------------------------------------------------


# lower bound at location parameter (false in lmom tutorial)

#' @title fitlm_GPD
#'
#' @description Function for fitting the Generalized Pareto Distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param bound Is the distribution bound? Default is NULL.
#'
#' @export
#'

fitlm_GPD=function(x,bound=NULL) {
  x <- na.omit(coredata(x))
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

  MLE<--sum(log(dgpd(x,location=fit$location,scale=fit$scale,shape=fit$shape)))
  plotpos<-lmomco::pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qgpd(p=plotpos,location=fit$location,scale=fit$scale,shape=fit$shape)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]

  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- list(CramerVonMises = CM,KolmogorovSmirnov = KS,MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'gpd')

  Res<-list()
  Res$Distribution<-list(FXs="qgpd")
  Res$Param<-fit
  Res$TheorLMom<-lmrgpa(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF = GoF

  return(Res)
}


# Generalised Gamma -------------------------------------------------------------

#' @name GenGamma
#' @title Generalized Gamma Distribution
#'
#' @description Generalized Gamma distribution
#'
#' @export
#'

dgengamma=function(x,scale, shape1, shape2){
  require(VGAM)
  fx = dgengamma.stacy(x = x,scale = scale, k = (shape1/shape2), d = shape2)
  return(fx)
}

#' @rdname GenGamma
#' @export
pgengamma=function(q,scale, shape1, shape2){
  require(VGAM)
  FX = pgengamma.stacy(q = q ,scale = scale, k = (shape1/shape2), d = shape2)
  return(FX)
}

#' @rdname GenGamma
#' @export
qgengamma=function(p,scale, shape1, shape2){
  require(VGAM)
  X = qgengamma.stacy(p =  p ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

#' @rdname GenGamma
#' @export
rgengamma=function(n,scale, shape1, shape2){
  require(VGAM)
  X = rgengamma.stacy(n =  n ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

#' @title fitlm_gengamma
#'
#' @description Function for fitting the Generalized Gamma distribution using the L-Moments method
#' Since there is not closed expression for the L-Moments of this distribution, the fitting must be done numerically.
#' For this purpose we employ the lmom::pelp function which works by optimization and requires an initial set of parameters.
#'
#' @param x A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param method The method to determine the starting values for the fitting. Can be 'knn' or 'DEoptim'. Default is 'knn'. See description for details.
#'
#' @export
#'

fitlm_gengamma=function(x,ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:5))  {
  max_order = max(4,order)
  x <- na.omit(coredata(x))
  PW = 1
  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }

  pfunction=anyFit::pgengamma
  qfunction=anyFit::qgengamma
  dfunction = anyFit::dgengamma

  sample_LM = Lmoments::Lcoefs(NZ, rmax = max_order, na.rm = FALSE, trim = c(0, 0))
  #I = RANN::nn2(Burr_InitValues[,c('L1','L2','tau3','tau4')],query = sample_LM,k = 10)$nn.idx
  init_values = GG_InitValues
  target_LMs = data.frame(lcv = sample_LM[2]/sample_LM[1], tau_3 = sample_LM[3], tau_4 = sample_LM[4])
  I = RANN::nn2(init_values[,c('lcv','tau_3','tau_4')],query = target_LMs,k = 50)$nn.idx
  start_matrix = init_values[I,]
  i=1
  params_optim <- function(params, target_LMs, order){
    max_order = max(order)
    scale = params[1]
    shape1 = params[2]
    shape2 = params[3]
    temp_lms <- lmrq(qfunction, order = c(1:max_order),
                     scale=scale, shape1=shape1, shape2=shape2,
                     subdiv = 10^6,acc = 10^-3)
    # temp_lms[5] =  temp_lms[2]/temp_lms[1]
    temp_err <- sapply(order, FUN = function(x){((target_LMs[x] - temp_lms[x])/target_LMs[x])^2})
    temp_err <- sqrt(sum(temp_err))
    return(temp_err)
  }
  all = optim(as.numeric(start_matrix[i,c(1,2,3)]), fn = params_optim, target_LMs = sample_LM, order = order,
              method = "L-BFGS-B", lower = c(1, 0.005, 0.05), upper = c(500,2000,200))
  params = list(scale = all$par[1], shape1 = all$par[2], shape2 = all$par[3])

  TheorLmom=lmrq(qfunction, order = c(1:5),
                 scale=params$scale, shape1=params$shape1, shape2=params$shape2,
                 subdiv = 10000,acc = 10^-2)

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
#' @description Generalized Gamma with Location distribution
#'
#' @export
#'

pgengamma_loc=function(q, location, scale, shape1, shape2){
  require(VGAM)
  FX = pgengamma.stacy(q = q-location, scale = scale, k = (shape1/shape2), d = shape2)
  return(FX)
}

#' @rdname GenGamma-Location
#' @export
dgengamma_loc=function(x,location, scale, shape1, shape2){
  require(VGAM)
  fx = dgengamma.stacy(x = x - location,scale = scale, k = (shape1/shape2), d = shape2)
  return(fx)
}

#' @rdname GenGamma-Location
#' @export
qgengamma_loc=function(p, location, scale, shape1, shape2){
  require(VGAM)
  X = location+ qgengamma.stacy(p =  p ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

#' @rdname GenGamma-Location
#' @export
rgengamma_loc=function(n, location, scale, shape1, shape2){
  require(VGAM)
  X = location + rgengamma.stacy(n =  n ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

#' @title fitlm_gengamma_loc
#'
#' @description Function for fitting the Generalized Gamma with location distribution using the L-Moments method.
#' Since there is not closed expression for the L-Moments of this distribution, the fitting must be done numerically.
#' For this purpose we employ the lmom::pelp function which works by optimization and requires an initial set of parameters.
#'
#' @param x A xts object containing the time series data.
#' @param location The location parameter of the distribution
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param method The method to determine the starting values for the fitting. Can be 'knn' or 'DEoptim'. Default is 'knn'. See description for details.
#'
#' @export
#'


fitlm_gengamma=function(x, location, ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:5))  {
  max_order = max(4,order)
  x <- na.omit(coredata(x))
  PW = 1
  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }

  sample_LM = Lmoments::Lcoefs(NZ, rmax = max_order, na.rm = FALSE, trim = c(0, 0))

  temp_fit = fitlm_gengamma(x = NZ - location, ignore_zeros = ignore_zeros, zero_threshold = zero_threshold, order = order)

  TheorLmom= temp_fit$TheorLmom
  TheorLmom$lambda_1 = TheorLmom$lambda_1 + location

  GoF <- GOF_tests(x = NZ, fit = temp_fit$params, distribution = 'gengamma_loc')

  #para$PW=PW

  Res<-list()
  Res$Distribution<-list(FXs="qgengamma_loc")
  Res$Param<-temp_fit$params
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
#' @description Burr Type XII distribution
#'
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
  p=runif(n)
  q=qburrDK2(p = p, scale=scale, shape1=shape1, shape2=shape2, PW=PW)
  return(q)
}

#' @title fitlm_burr
#'
#' @description Function for fitting the Burr Type XII distribution using the L-Moments method
#' Since there is not closed expression for the L-Moments of this distribution, the fitting must be done numerically.
#' For this purpose we employ the lmom::pelp function which works by optimization and requires an initial set of parameters.
#'
#' @param x A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param method The method to determine the starting values for the fitting. Can be 'knn' or 'DEoptim'. Default is 'knn'. See description for details.
#'
#' @export
#'

fitlm_burr=function(x,ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:5))  {
  max_order = max(4,order)
  x <- na.omit(coredata(x))
  PW = 1
  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }

  pfunction=anyFit::pburr
  qfunction=anyFit::qburr
  dfunction = anyFit::dburr

  sample_LM = Lmoments::Lcoefs(NZ, rmax = max_order, na.rm = FALSE, trim = c(0, 0))
  #I = RANN::nn2(Burr_InitValues[,c('L1','L2','tau3','tau4')],query = sample_LM,k = 10)$nn.idx
  init_values = Burr_InitValues
  target_LMs = data.frame(lcv = sample_LM[2]/sample_LM[1], tau_3 = sample_LM[3], tau_4 = sample_LM[4])
  I = RANN::nn2(init_values[,c('lcv','tau_3','tau_4')],query = target_LMs,k = 50)$nn.idx
  start_matrix = init_values[I,]
  i=1
  params_optim <- function(params, target_LMs, order){
    max_order = max(order)
    scale = params[1]
    shape1 = params[2]
    shape2 = params[3]
    temp_lms <- lmrp(pfunction, bounds = c(0,Inf),order = c(1:max_order),
                     scale=scale, shape1=shape1, shape2=shape2,
                     subdiv = 10000,acc = 10^-3)
    # temp_lms[5] =  temp_lms[2]/temp_lms[1]
    temp_err <- sapply(order, FUN = function(x){((target_LMs[x] - temp_lms[x])/target_LMs[x])^2})
    temp_err <- sqrt(sum(temp_err))
    return(temp_err)
  }
  all = optim(as.numeric(start_matrix[i,c(1,2,3)]), fn = params_optim, target_LMs = sample_LM, order = order,
              method = "L-BFGS-B", lower = c(0.5, 0.5, 0.001), upper = c(50,50,1))
  params = list(scale = all$par[1], shape1 = all$par[2], shape2 = all$par[3])

  TheorLmom=lmrp(pfunction, bounds = c(0, Inf), order = c(1:5),
                 scale=params$scale, shape1=params$shape1, shape2=params$shape2,
                 subdiv = 10000,acc = 10^-2)

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
#' @description Dagum distribution
#'
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
  qdagum(p = runif(n), scale=scale, shape1 = shape1, shape2 = shape2, PW=PW)
}


#' @title fitlm_dagum
#'
#' @description Function for fitting the Dagum distribution using the L-Moments method.
#' Since there is not closed expression for the L-Moments of this distribution, the fitting must be done numerically.
#' For this purpose we employ the lmom::pelp function which works by optimization and requires an initial set of parameters.
#'
#' @param x A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param method The method to determine the starting values for the fitting. Can be 'knn' or 'DEoptim'. Default is 'knn'. See description for details.
#'
#' @export
#'

fitlm_dagum=function(x,ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:5))  {
  max_order = max(4,order)
  x <- na.omit(coredata(x))
  PW = 1
  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }

  pfunction=anyFit::pdagum
  qfunction=anyFit::qdagum
  dfunction = anyFit::ddagum

  sample_LM = Lmoments::Lcoefs(NZ, rmax = max_order, na.rm = FALSE, trim = c(0, 0))
  #I = RANN::nn2(Burr_InitValues[,c('L1','L2','tau3','tau4')],query = sample_LM,k = 10)$nn.idx
  init_values = Dagum_InitValues
  target_LMs = data.frame(lcv = sample_LM[2]/sample_LM[1], tau_3 = sample_LM[3], tau_4 = sample_LM[4])
  I = RANN::nn2(init_values[,c('lcv','tau_3','tau_4')],query = target_LMs,k = 50)$nn.idx
  start_matrix = init_values[I,]
  i=1
  params_optim <- function(params, target_LMs, order){
    max_order = max(order)
    scale = params[1]
    shape1 = params[2]
    shape2 = params[3]
    temp_lms <- lmrp(pfunction, bounds = c(0,Inf),order = c(1:max_order),
                     scale=scale, shape1=shape1, shape2=shape2,
                     subdiv = 10000,acc = 10^-3)
    # temp_lms[5] =  temp_lms[2]/temp_lms[1]
    temp_err <- sapply(order, FUN = function(x){((target_LMs[x] - temp_lms[x])/target_LMs[x])^2})
    temp_err <- sqrt(sum(temp_err))
    return(temp_err)
  }
  all = optim(as.numeric(start_matrix[i,c(1,2,3)]), fn = params_optim, target_LMs = sample_LM, order = order,
              method = "L-BFGS-B", lower = c(0.5, 0.0001, 0.000001), upper = c(1000,500,0.5))
  params = list(scale = all$par[1], shape1 = all$par[2], shape2 = all$par[3])

  TheorLmom=lmrp(pfunction, bounds = c(0, Inf), order = c(1:5),
                 scale=params$scale, shape1=params$shape1, shape2=params$shape2,
                 subdiv = 10000,acc = 10^-2)

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
#' @description Exponential Weibull distribution
#'
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
  u = runif(n)
  sim <-  stats::qweibull(u^(1/shape2),scale=scale,shape=shape1)
  return(sim)
}

#' @title fitlm_expweibull
#'
#' @description Function for fitting the Exponential Weibull distribution using the L-Moments method
#'
#' @param x A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @export
#'

fitlm_expweibull=function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  # Important: MEANT to be fitted ONLY to NON-NEGATIVE data.
  # The potential inlcusion of zeros or NAs to "x" is hanlded by the 2 lines below.
  xNZ <- na.omit(coredata(x))

  if (ignore_zeros == TRUE){
    xNZ <- x[x > zero_threshold,]
    xNZ <- na.omit(coredata(xNZ))
  }
  lm=lmom::samlmu(xNZ)

  # IDEA for helping pelq with "good" statring values
  # Step1: Fit the classic Weibull distribution (scale, shape)
  # Step2: Use the obtained (scale, shape) as starting values for the pelq function.
  #       In detail, set:
  #                      scale=scale (as found from fitting the classic Weibull in Step1)
  #                      shape1=shape (as found from fitting the classic Weibull in Step1)
  #                      shape2=1 (start from 1, since it correspond to the classic Weibull)
  #
  # The above strategy is working fine in most cases, yet some fine-tuning may be needed to shape2
  # (e.g., using try) in case of pelq returns an error.

  # Implement Step1
  parwei=lmom::pelwei(lmom = lm, bound = 0)[2:3]

  # Implement Step2
  START=c(parwei[1], parwei[2], 1)
  fit=lmom::pelq(lmom = lm[1:3], qfunc = qexpweibull, start = START, type = 's')$para
  fit=list(scale=as.vector(fit[1]), shape1=as.vector(fit[2]), shape2=as.vector(fit[3]))

  GoF <- GOF_tests(x = x, fit = fit, distribution = 'expweibull')

  Res<-list()
  Res$Distribution<-list(FXs="qexpweibull")
  Res$Param<-fit
  Res$GoF<-GoF

  return(Res)

}

