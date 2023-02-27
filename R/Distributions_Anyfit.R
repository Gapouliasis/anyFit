#' @title anyFit_distributions
#' 
#' @description anyFit distributions and distribution fitting functions 
#' 



#Anyfit supported distributions
library("lmom")
library("lmomco")

library("fitdistrplus")
library("PearsonDS")
library("FAdist")
library("EnvStats")
library("VGAM")
library("ExtDist")

# Exponential -------------------------------------------------------------

# x<-rexp(1000,location=0,scale=3) # my function

dexp=function(x,location,scale) {
  fx=stats::dexp(x=x-location,rate=1/scale)
  return(fx)
}

pexp=function(q,location,scale) {
  FX=lmom::cdfexp(x=q,para=c(location,scale))
  return(FX)
}

qexp=function(p,location,scale) {
  x=lmom::quaexp(f=p,para=c(location,scale))
  return(x)
}

rexp=function(n,location,scale) {
  x=lmom::quaexp(f=runif(n),para=c(location,scale))
  return(x)
}


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
  u_emp <- ppoints(x)
  q_emp <- x
  qq_fitted <- qexp(u_emp, location = fit$location, scale = fit$scale)
  
  CM <- CDFt::CramerVonMisesTwoSamples(q_emp, qq_fitted)
  KS <- CDFt::KolmogorovSmirnov(q_emp, qq_fitted)
  
  MLE<--sum(log(dexp(x,location=fit$location,scale=fit$scale)))
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qexp(p=plotpos,location=fit$location,scale=fit$scale)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  
  GoF <- list(MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
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

drayleigh = function(x, location, scale) {
  fx <- VGAM::drayleigh(x - location, scale, log = FALSE)
  return(fx)
}

prayleigh = function(x, location, scale) {
  fx <- VGAM::prayleigh(x - location, scale, lower.tail = TRUE, log.p = FALSE)  
  return(fx)
}

qrayleigh = function(p, location, scale) {
  fx <- VGAM::qrayleigh(p, scale, lower.tail = TRUE, log.p = FALSE) + location
  return(fx)
}

rrayleigh = function(n, location, scale) {
  fx <- VGAM::rrayleigh(n, scale) + location
  return(fx)
}

fitlm_rayleigh = function(x,ignore_zeros = FALSE, zero_threshold = 0.01) {
  x <- na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }
  l1 <- lmom::samlmu(x,nmom = 1, ratios = FALSE)
  l2 <- lmom::samlmu(x,nmom = 2, ratios = FALSE)
  location <- l1 - sqrt(2)/(sqrt(2) - 1)*l2
  ratio <- (gamma(1.5))^2*(3-2*sqrt(2))/(2*l2^2)
  scale <- 1/sqrt(2*ratio)
  par <- c(location,scale)
  
  u_emp <- ppoints(x)
  q_emp <- x
  qq_fitted <- qrayleigh(u_emp, location = fit$location, scale = fit$scale)
  
  CM <- CDFt::CramerVonMisesTwoSamples(q_emp, qq_fitted)
  KS <- CDFt::KolmogorovSmirnov(q_emp, qq_fitted)
  MLE<--sum(log(dexp(x,location=location,scale=scale)))
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qrayleigh(p=plotpos,location=location,scale=scale)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  
  GoF <- list(CramerVonMises = CM,KolmogorovSmirnov = KS,MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  Res<-list()
  Res$Distribution<-list(FXs="qrayleigh")
  Res$Param<-par
  Res$TheorLMom<-lmrexp(par,nmom=4)
  names(sam)<-c('lambda_1','lambda_2','tau_3','tau_4')      
  Res$DataLMom<-sam
  Res$GoF<-GoF
  
  return(Res)
}



# Gamma -------------------------------------------------------------

# x<-stats::rgamma(10000,scale=3,shape=0.5)

dgamma=function(x,scale,shape) {
  fx=stats::dgamma(x=x,shape=shape,scale=scale)
  return(fx)
}

pgamma=function(q,scale,shape) {
  FX=stats::pgamma(q=q,shape=shape,scale=scale)
  return(FX)
}

qgamma=function(p,scale,shape) {
  x=stats::qgamma(p=p,shape=shape,scale=scale)
  return(x)
}

rgamma=function(n,scale,shape) {
  x=stats::rgamma(n=n,shape=shape,scale=scale)
  return(x)
}

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
  
  
  MLE<--sum(log(dgamma(x,scale=fit$scale,shape=fit$shape)))
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qgamma(p=plotpos,scale=fit$scale,shape=fit$shape)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- c(MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
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

dgamma3=function(x,location,scale,shape) {
  fx=PearsonDS::dpearsonIII(x=x,shape=shape,location=location,scale=1/scale)
  return(fx)
}

pgamma3=function(q,location,scale,shape) {
  FX=PearsonDS::ppearsonIII(q=q,shape=shape,location=location,scale=1/scale)
  return(FX)
}

qgamma3=function(p,location,scale,shape) {
  x=PearsonDS::qpearsonIII(p=p,shape=shape,location=location,scale=1/scale)
  return(x)
}

rgamma3=function(n,location,scale,shape) {
  x=PearsonDS::rpearsonIII(n=n,shape=shape,location=location,scale=1/scale)
  return(x)
}

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
  
  
  MLE<--sum(log(dgamma3(x,location=fit$location,scale=fit$scale,shape=fit$shape)))
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qgamma3(p=plotpos,location=fit$location,scale=fit$scale,shape=fit$shape)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- list(MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
  Res<-list()
  Res$Distribution<-list(FXs="qgamma3")
  Res$Param<-fit
  Res$TheorLMom<-lmrpe3(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF<-GoF
  
  return(Res)
  
}


# Generalised Logistic -------------------------------------------------------------

# x<-quaglo(runif(10000), c(0,1,-0.5))

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
  
  GoF <- list(CramerVonMises = CM,KolmogorovSmirnov = KS,MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
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

dnorm=function(x,mean=0,sd=1) {
  fx=stats::dnorm(x=x,sd=sd,mean=mean)
  return(fx)
}

pnorm=function(q,mean=0,sd=1) {
  FX=stats::pnorm(q=q,sd=sd,mean=mean)
  return(FX)
}

qnorm=function(p,mean=0,sd=1) {
  x=stats::qnorm(p=p,sd=sd,mean=mean)
  return(x)
}

rnorm=function(n,mean=0,sd=1) {
  x=stats::rnorm(n=n,sd=sd,mean=mean)
  return(x)
}


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
  
  MLE<--sum(log(dnorm(x,mean=fit$mean,sd=fit$sd)))
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qnorm(p=plotpos,mean=fit$mean,sd=fit$sd)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- list(CramerVonMises = CM,KolmogorovSmirnov = KS,MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
  Res<-list()
  Res$Distribution<-list(FXs="qnorm")
  Res$Param<-fit
  Res$TheorLMom<-lmrnor(par,nmom=4)
  Res$DataLMom<-sam
  Res$GoF<-GoF
  
  return(Res)
}



# Weibull, 3-parameter -------------------------------------------------------------

# x<-quawei(runif(10000), c(0,5,1.5)) # require(lmom)
# x<-rweibull(10,scale = 5,shape = 1.5)
# identical runs

dweibull=function(x,location,scale,shape) {
  fx=FAdist::dweibull3(x=x,shape=shape,scale=scale,thres=location)
  return(fx)
}

pweibull=function(q,location,scale,shape) {
  FX=FAdist::pweibull3(q=q,shape=shape,scale=scale,thres=location)
  return(FX)
}

qweibull=function(p,location,scale,shape) {
  x=FAdist::qweibull3(p=p,shape=shape,scale=scale,thres=location)
  return(x)
}

rweibull=function(n,location,scale,shape) {
  x=FAdist::rweibull3(n=n,shape=shape,scale=scale,thres=location)
  return(x)
}

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
  
  
  # print(paste("Weibull distribution bounded at",fit$location,sep=" "))
  
  MLE<--sum(log(dweibull(x,location=fit$location,scale=fit$scale,shape=fit$shape)))
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qweibull(p=plotpos,location=fit$location,scale=fit$scale,shape=fit$shape)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- list(MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
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

dgumbel=function(x,location,scale) {
  fx=FAdist::dgumbel(x=x,scale=scale,location=location)
  return(fx)
}

pgumbel=function(q,location,scale) {
  FX=FAdist::pgumbel(q=q,scale=scale,location=location)
  return(FX)
}

qgumbel=function(p,location,scale) {
  x=FAdist::qgumbel(p=p,scale=scale,location=location)
  return(x)
}

rgumbel=function(n,location,scale) {
  x=FAdist::rgumbel(n=n,scale=scale,location=location)
  return(x)
}


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
  
  
  MLE<--sum(log(dgumbel(x,location=fit$location,scale=fit$scale)))
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qgumbel(p=plotpos,location=fit$location,scale=fit$scale)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- list(CramerVonMises = CM,KolmogorovSmirnov = KS,MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
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

dlognorm=function(x,location=0,scale,shape) {
  fx=FAdist::dlnorm3(x=x,shape=shape,scale=scale,thres=location)
  return(fx)
}

plognorm=function(q,location=0,scale,shape) {
  FX=FAdist::plnorm3(q=q,shape=shape,scale=scale,thres=location)
  return(FX)
}

qlognorm=function(p,location=0,scale,shape) {
  x=FAdist::qlnorm3(p=p,shape=shape,scale=scale,thres=location)
  return(x)
}

rlognorm=function(n,location=0,scale,shape) {
  x=FAdist::rlnorm3(n=n,shape=shape,scale=scale,thres=location)
  return(x)
}

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
  
  
  # print(paste("Lognormal distribution bounded at",fit$location,sep=" "))
  
  MLE<--sum(log(dlnorm(x,location=fit$location,scale=fit$scale,shape=fit$shape)))
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qlnorm(p=plotpos,location=fit$location,scale=fit$scale,shape=fit$shape)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- list(CramerVonMises = CM,KolmogorovSmirnov = KS,MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
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

dgev=function(x,location,scale,shape) {
  fx=FAdist::dgev(x=x,shape=shape,scale=scale,location=location)
  return(fx)
}

pgev=function(q,location,scale,shape) {
  FX=FAdist::pgev(q=q,shape=shape,scale=scale,location=location)
  return(FX)
}

qgev=function(p,location,scale,shape) {
  x=FAdist::qgev(p=p,shape=shape,scale=scale,location=location)
  return(x)
}

rgev=function(n,location,scale,shape) {
  x=FAdist::rgev(n=n,shape=shape,scale=scale,location=location)
  return(x)
}


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
  
  MLE<--sum(log(dgev(x,location=fit$location,scale=fit$scale,shape=fit$shape)))
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qgev(p=plotpos,location=fit$location,scale=fit$scale,shape=fit$shape)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- list(CramerVonMises = CM,KolmogorovSmirnov = KS,MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
  Res<-list()
  Res$Distribution<-list(FXs="qgev")
  Res$Param<-fit
  Res$TheorLMom<-lmrgev(par,nmom=4)
  Res$DataLMom<-sam
  
  return(Res)
}



# GenPareto (need change) -------------------------------------------------------------


# lower bound at location parameter (false in lmom tutorial)

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
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qgpd(p=plotpos,location=fit$location,scale=fit$scale,shape=fit$shape)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- list(CramerVonMises = CM,KolmogorovSmirnov = KS,MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
  
  Res<-list()
  Res$Param<-fit
  Res$TheorLMom<-lmrgpa(par,nmom=4)
  Res$DataLMom<-sam
  
  return(Res)
}


# Generalised Gamma -------------------------------------------------------------

dgengamma=function(x,scale, shape1, shape2){
  require(VGAM)
  fx = dgengamma.stacy(x = x,scale = scale, k = (shape1/shape2), d = shape2)
  return(fx)
}

pgengamma=function(q,scale, shape1, shape2){
  require(VGAM)
  FX = pgengamma.stacy(q = q ,scale = scale, k = (shape1/shape2), d = shape2)
  return(FX)
}

qgengamma=function(p,scale, shape1, shape2){
  require(VGAM)
  X = qgengamma.stacy(p =  p ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

rgengamma=function(n,scale, shape1, shape2){
  require(VGAM)
  X = rgengamma.stacy(n =  n ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

fitlm_gengamma=function(x,ignore_zeros = FALSE, zero_threshold = 0.01)  {
  x <- na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }
  
  lm=lmom::samlmu(NZ)[1:3]
  pfunction=pgengamma
  qfunction=qgengamma
  dfunction = dgengamma
  MLE_fun = function(trial = c(par1,par2,par3),x_ts,dfunction){
    #trial = c(par1,par2,par3)
    dist_args = formalArgs(dfunction)
    temp_args = list(x_ts)
    temp_args = c(temp_args, as.list(trial))
    names(temp_args) = dist_args
    LF = do.call(dfunction, temp_args)
    LF = -sum(log(LF))
    if (is.infinite(LF)){LF = 10^10}
    return(LF)
  }
  
  itmax = 0
  NP = 10*(length(formalArgs(dfunction)) - 1)
  for (i in 1:10){
    test_value <- tryCatch({
      itmax = 40
      start = DEoptim::DEoptim(MLE_fun, lower = c(0.001,0.001,0.001), upper = c(30,20,5), x_ts= NZ, 
                               dfunction = dfunction, control=list(itermax=itmax, NP = NP,trace = FALSE,
                                                                   F = 0.65, parallelType=1))
      start = start$optim$bestmem
      para = pelp(lmom = lm,
                  pfunc = pfunction,
                  start = start,
                  bounds = c(0, Inf),
                  type = 's')
      para = as.list(para$para)},
      warning=function(cond) {
        return(para)},
      error = function(cond) { 
        return(NA)
        message(cond)})
    #}
    
    
    if (is.na(test_value[1])){
      itmax = itmax + 10
      NP = NP + 10
    }else{
      TheorLmom=lmrp(pfunction, bounds = c(0, Inf), order = c(1:3),
                     scale=para$scale, shape1=para$shape1, shape2=para$shape2,
                     subdiv = 10000,acc = 10^-2)
      if (abs(lm[1]-TheorLmom[1])>1000){
        itmax = itmax + 10
        NP = NP + 10
      }else{
        para = test_value
        break
      }
    }
  }
  

  TheorLmom=lmrp(pfunction, bounds = c(0, Inf), order = c(1:3),
                 scale=para$scale, shape1=para$shape1, shape2=para$shape2,
                 subdiv = 10000,acc = 10^-2)

  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qgengamma(p=plotpos,scale=para$scale, shape1=para$shape1, shape2=para$shape2)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- c(MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
  
  #para$PW=PW
  
  Res<-list()
  Res$Distribution<-list(FXs="qgengamma")
  Res$Param<-para
  Res$TheorLMom<-TheorLmom
  Res$DataLMom<-lm
  Res$GoF<-GoF
  Res$Diag = c(itmax,NP)
  
  return(Res)
}


#Generalized Gamma with location----------------------------------------

pgengamma_loc=function(q, location, scale, shape1, shape2){
  require(VGAM)
  FX = pgengamma.stacy(q = q-location, scale = scale, k = (shape1/shape2), d = shape2)
  return(FX)
}

dgengamma_loc=function(x,location, scale, shape1, shape2){
  require(VGAM)
  fx = dgengamma.stacy(x = x - location,scale = scale, k = (shape1/shape2), d = shape2)
  return(fx)
}

qgengamma_loc=function(p, location, scale, shape1, shape2){
  require(VGAM)
  X = location+ qgengamma.stacy(p =  p ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

rgengamma_loc=function(n, location, scale, shape1, shape2){
  require(VGAM)
  X = location + rgengamma.stacy(n =  n ,scale = scale, k = (shape1/shape2), d = shape2)
  return(X)
}

fitlm_gengamma_loc=function(x,location,ignore_zeros = FALSE, zero_threshold = 0.01)  {
  x <- na.omit(coredata(x))
  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }
  
  lm=lmom::samlmu(NZ)[1:3]
  pfunction=pgengamma_loc
  qfunction=qgengamma_loc
  dfunction = dgengamma_loc
  MLE_fun = function(trial = c(par1,par2,par3),x_ts,dfunction,location){
    #trial = c(par1,par2,par3)
    dist_args = formalArgs(dfunction)
    temp_args = list(x_ts)
    temp_args = c(temp_args, as.list(trial))
    names(temp_args) = dist_args
    LF = do.call(dfunction, temp_args)
    LF = -sum(log(LF))
    if (is.infinite(LF)){LF = 10^10}
    return(LF)
  }
  
  itmax = 0
  NP = 10*(length(formalArgs(dfunction)) - 1)
  for (i in 1:10){
    test_value <- tryCatch({
      itmax = 40
      start = DEoptim::DEoptim(MLE_fun, lower = c(0.001,0.001,0.001), upper = c(30,20,5), x_ts= NZ, 
                               dfunction = dfunction, control=list(itermax=itmax, NP = NP,trace = FALSE,
                                                                   F = 0.65, parallelType=1))
      start = start$optim$bestmem
      para = pelp(lmom = lm,
                  pfunc = pfunction,
                  start = start,
                  bounds = c(0, Inf),
                  type = 's')
      para = as.list(para$para)},
      warning=function(cond) {
        return(para)},
      error = function(cond) { 
        return(NA)
        message(cond)})
    #}
    
    
    if (is.na(test_value[1])){
      itmax = itmax + 10
      NP = NP + 10
    }else{
      TheorLmom=lmrp(pfunction, bounds = c(0, Inf), order = c(1:3),
                     scale=para$scale, shape1=para$shape1, shape2=para$shape2,
                     subdiv = 10000,acc = 10^-2)
      if (abs(lm[1]-TheorLmom[1])>1000){
        itmax = itmax + 10
        NP = NP + 10
      }else{
        para = test_value
        break
      }
    }
  }
  
  
  TheorLmom=lmrp(pfunction, bounds = c(0, Inf), order = c(1:3),
                 scale=para$scale, shape1=para$shape1, shape2=para$shape2,
                 subdiv = 10000,acc = 10^-2)
  
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qgengamma(p=plotpos,scale=para$scale, shape1=para$shape1, shape2=para$shape2)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- c(MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
  
  #para$PW=PW
  
  Res<-list()
  Res$Distribution<-list(FXs="qgengamma")
  Res$Param<-para
  Res$TheorLMom<-TheorLmom
  Res$DataLMom<-lm
  Res$GoF<-GoF
  Res$Diag = c(itmax,NP)
  
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
dburr=function(q, scale, shape1, shape2, PW=1){
  d=PW*shape1*scale^(-shape1)*q^(shape1-1)*((shape1*shape2*(q/scale)^shape1)+1)^(-1/(shape1*shape2)-1)
}

pburr=function(q, scale, shape1, shape2, PW=1) {
  # scale: lamda
  # shape1: zeta
  # shape2: tail index!
  p=1-PW*(1+(shape1*shape2)*(q/scale)^shape1)^(-1/(shape1*shape2))
  return(p)
}

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

rburr=function(n, scale, shape1, shape2, PW=1) {
  p=runif(n)
  q=qburrDK2(p = p, scale=scale, shape1=shape1, shape2=shape2, PW=PW)
  return(q)
}

fitlm_burr=function(x,ignore_zeros = FALSE, zero_threshold = 0.01)  {
  x <- na.omit(coredata(x))
  PW = 1
  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }
  
  lm=lmom::samlmu(NZ)[1:3]
  pfunction=pburr
  qfunction=qburr
  dfunction = dburr
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
  
  itmax = 0
  NP = 10*(length(formalArgs(dfunction)) - 1)
  for (i in 1:10){
    test_value <- tryCatch({
      itmax = 40
      start = DEoptim::DEoptim(MLE_fun, lower = c(0.001,0.001,0.001), upper = c(40,40,40), x_ts= NZ, 
                               dfunction = dfunction, control=list(itermax=itmax, NP = NP,trace = FALSE,
                                                                   F = 0.65, parallelType=1))
      start = start$optim$bestmem
      para = pelp(lmom = lm,
                  pfunc = pfunction,
                  start = start,
                  bounds = c(0, Inf),
                  type = 's')
      para = as.list(para$para)},
      warning=function(cond) {
        return(para)},
      error = function(cond) { 
        return(NA)
        message(cond)})
    if (is.na(test_value[1])){
      itmax = itmax + 10
      NP = NP + 10
    }else{
      TheorLmom=lmrp(pfunction, bounds = c(0, Inf), order = c(1:3),
                     scale=para$scale, shape1=para$shape1, shape2=para$shape2,
                     subdiv = 10000,acc = 10^-2)
      if (abs((lm[1]-TheorLmom[1])/lm[1])>0.1){
        itmax = itmax + 10
        NP = NP + 10
      }else{
        para = test_value
        break
      }
    }
  }
  
  
  TheorLmom=lmrp(pfunction, bounds = c(0, Inf), order = c(1:3),
                 scale=para$scale, shape1=para$shape1, shape2=para$shape2,
                 subdiv = 10000,acc = 10^-2)
  
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qburr(p=plotpos,scale=para$scale, shape1=para$shape1, shape2=para$shape2)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- c(MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
  
  #para$PW=PW
  
  Res<-list()
  Res$Distribution<-list(FXs="qburr")
  Res$Param<-para
  Res$TheorLMom<-TheorLmom
  Res$DataLMom<-lm
  Res$GoF<-GoF
  Res$Diag = c(itmax,NP)
  
  return(Res)
}


#Dagum--------------------------------------------------------------
pdagum=function(q, scale, shape1, shape2, PW=1) {
  # The moments exist (i.e., have finite values) only for order r < 1/shape2;  for larger r they are infinite.
  p=(1-PW)+PW*(1+(shape2/shape1)*(q/scale)^(-1/shape2))^-shape1
  return(p)
}

ddagum=function(q, scale, shape1, shape2, PW=1) {
  d=PW*(1/scale)*(q/scale)^(-1-1/shape2)*((1/shape1)*(shape1 + shape2*(q/scale)^(-1/shape2)))^(-1-shape1)
}

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

rdagum=function(n, scale, shape1, shape2, PW=1) {
  qdagum(p = runif(n), scale=scale, shape1 = shape1, shape2 = shape2, PW=PW)
}

fitlm_dagum=function(x,ignore_zeros = FALSE, zero_threshold = 0.01)  {
  x <- na.omit(coredata(x))
  PW = 1
  if (ignore_zeros == TRUE){
    NZ=x[x>zero_threshold,]
  }else{
    NZ=x
  }
  
  lm=lmom::samlmu(NZ)[1:3]
  pfunction=pdagum
  qfunction=qdagum
  dfunction = ddagum
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
  
  itmax = 0
  NP = 10*(length(formalArgs(dfunction)) - 1)
  for (i in 1:10){
    test_value <- tryCatch({
      itmax = 40
      start = DEoptim::DEoptim(MLE_fun, lower = c(0.001,0.001,0.001), upper = c(40,40,40), x_ts= NZ, 
                               dfunction = dfunction, control=list(itermax=itmax, NP = NP,trace = FALSE,
                                                                   F = 0.65, parallelType=1))
      start = start$optim$bestmem
      para = pelp(lmom = lm,
                  pfunc = pfunction,
                  start = start,
                  bounds = c(0, Inf),
                  type = 's')
      para = as.list(para$para)},
      warning=function(cond) {
        return(para)},
      error = function(cond) { 
        return(NA)
        message(cond)})
    if (is.na(test_value[1])){
      itmax = itmax + 10
      NP = NP + 10
    }else{
      TheorLmom=lmrp(pfunction, bounds = c(0, Inf), order = c(1:3),
                     scale=para$scale, shape1=para$shape1, shape2=para$shape2,
                     subdiv = 10000,acc = 10^-2)
      if (abs((lm[1]-TheorLmom[1])/lm[1])>0.1){
        itmax = itmax + 10
        NP = NP + 10
      }else{
        para = test_value
        break
      }
    }
  }
  
  TheorLmom=lmrp(pfunction, bounds = c(0, Inf), order = c(1:3),
                 scale=para$scale, shape1=para$shape1, shape2=para$shape2,
                 subdiv = 10000,acc = 10^-2)
  
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qdagum(p=plotpos,scale=para$scale, shape1=para$shape1, shape2=para$shape2)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- c(MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
  
  #para$PW=PW
  
  Res<-list()
  Res$Distribution<-list(FXs="qdagum")
  Res$Param<-para
  Res$TheorLMom<-TheorLmom
  Res$DataLMom<-lm
  Res$GoF<-GoF
  Res$Diag = c(itmax,NP)
  
  return(Res)
}

# Exponentiated Weibull (Generalized Weibull) Distribution-----------------------------
# Cumulative distribution function
pexpweibull<- function(q,scale,shape1,shape2,log.p=FALSE){
  log.cdf <- shape2*stats::pweibull(q = q, scale=scale,shape=shape1,log.p=TRUE)
  ifelse(log.p, return(log.cdf), return(exp(log.cdf)))
}  

# Probability density function
dexpweibull<- function(x,scale,shape1,shape2,log=FALSE){
  log.pdf <-  log(shape2) + (shape2-1)*stats::pweibull(q = x,scale=scale,shape=shape1,log=TRUE) + 
    stats::dweibull(x = x,scale=scale,shape=shape1,log=TRUE)
  ifelse(log, return(log.pdf), return(exp(log.pdf)))
}

# Quantile function
qexpweibull<- function(p,scale,shape1,shape2){
  quant <-  stats::qweibull(p = p^(1/shape2),scale=scale,shape=shape1)
  return(quant)
}  

# Random number generation
rexpweibull<- function(n,scale,shape1,shape2){
  u = runif(n)
  sim <-  stats::qweibull(u^(1/shape2),scale=scale,shape=shape1)
  return(sim)
} 



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
  
  MLE<--sum(log(dexpweibull(x=xNZ,scale=fit$scale,shape1=fit$shape1,shape2=fit$shape2)))
  plotpos<-pp(x=x,a=0,sort=FALSE)
  theorquantiles<-qexpweibull(p=plotpos,scale=fit$scale,shape1=fit$shape1,shape2=fit$shape2)
  theorquantilessort<-sort(theorquantiles,decreasing=TRUE)[1:10]
  samplesort<-sort(x,decreasing=TRUE)[1:10]
  
  MSEquant<-sum((theorquantiles-x)^2)/length(x)
  DiffOfMax<-((max(theorquantiles)-max(x))/max(x))*100
  MeanDiffOf10Max<-sum(abs(theorquantilessort-samplesort))/10
  GoF <- list(MLE=MLE,
              MSEquant=MSEquant,DiffOfMax=DiffOfMax,MeanDiffOf10Max=MeanDiffOf10Max)
  
  Res<-list()
  Res$Distribution<-list(FXs="qexpweibull")
  Res$Param<-fit
  Res$GoF<-GoF
  
  return(Res)
  
}

