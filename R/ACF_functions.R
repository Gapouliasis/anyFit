#' @title ACF_functions
#'
#' @description  Theoretical ACF functions. 

CAS_ACF <- function(kappa,beta,lag_max = 10){
  p_cas <- data.frame(lag = seq(0,lag_max))
  p_cas$ACF <- (1+kappa*beta*p_cas$lag)^(-1/beta)
  return(p_cas)
}

HK_ACF <- function(H , lag_max = 10){
  p_hk = data.frame(lag = seq(0,lag_max))
  p_hk$ACF <- 0.5*((abs(p_hk$lag-1))^(2*H)-
                     2*(abs(p_hk$lag))^(2*H)+
                     (abs(p_hk$lag+1))^(2*H))
  return(p_hk)
}

SRD_ACF <- function(kappa, lag_max = 10){
  psrd = data.frame(lag = seq(0,lag_max))
  psrd$ACF <- exp(-psrd$lag*kappa)
  return(psrd)
}
