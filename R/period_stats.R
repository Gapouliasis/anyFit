#' @title period_stats
#' 
#' @description Function to calculate period based basic statistics. Specifically, the number of data points, 
#' number of missing data, percentage of missing data, min, max, mean, variance, 
#' standard deviation, variation, 3rd moment, skewness, kurtosis, l-statistics, i.e. mean, scale, 
#' 3rd and 4th order l-moments, l-moment variation, l-moment skewness, l-moment kurtosis, 
#' quantiles of 5, 25, 50, 75 and 95, and the inter quartile range. 
#' 
#'
#' @param ts a xts object containing the time series data 
#' @param period a period for which to calculate statistics. Default is NA
#' @param period_multiplier a multiplier to define arbitrary periods based on the period argument. 
#' E.g. 2 months, 3 years etc. Default is 1.
#'
#' @return a xts object containing the period-based statistics for the time series data. 
#' 
#' @examples
#'file <- "KNMI_Daily.csv"
#'file_path <- file.path(cwd,file)
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path, 
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'                  
#' period_stats(ts, period = "months", period_multiplier = 3)
#'
#' @export
#'

period_stats <-function(ts, period = NA, period_multiplier = 1){
  spec_period <- endpoints(ts, on = 'months', k = period_multiplier)
  
  NumofData <- period.apply(ts, spec_period, FUN = length)
  NA_vals <- xts(x = is.na(ts), order.by = index(ts))
  NumofMisData <- period.apply(NA_vals, spec_period, FUN = sum)
  PercOfMissingData <- (NumofMisData/NumofData)*100
  Min <- round(period.apply(ts, spec_period, FUN = min),2)
  Max <- round(period.apply(ts, spec_period, FUN = max),2)
  Mean <- round(period.apply(ts, spec_period, FUN = mean),2)
  Var <- round(period.apply(ts, spec_period, FUN = var),2)
  StDev <- round(period.apply(ts, spec_period, FUN = sd),2)
  Min <- round(period.apply(ts, spec_period, FUN = min),2)
  Variation <- round(StDev/Mean,2)
  Mom3 <- round(sum((ts-Mean)^3,na.rm=T)*(1/(sum(!is.na(ts))-1)),2)
  Skewness <- round(period.apply(ts, spec_period, FUN = moments::skewness),2)
  Kurtosis <- round(period.apply(ts, spec_period, FUN = moments::kurtosis),2)
  
  #L-moments 
  lmom <- round(period.apply(ts, spec_period, FUN = lmom::samlmu, nmom = 4, ratios = FALSE, trim = 0),2)
  LMean<-round(lmom[,1],2)
  LScale<-round(lmom[,2],2)
  L3<-round(lmom[,3],2)
  L4<-round(lmom[,4],2)
  lmom_ratios<-round(period.apply(ts, spec_period, FUN = lmom::samlmu, nmom = 4, ratios = TRUE, trim = 0),2)
  LVariation<-round(LScale/LMean,2)
  LSkewness<-round(lmom_ratios[,3],2)
  Lkurtosis<-round(lmom_ratios[,4],2)
  #Pdr <- round(mean(ts<=0.01,na.rm=T),2)
  
  # Quantiles of all values
  Q5 <- round(period.apply(ts, spec_period, FUN = quantile,probs=c(0.05),na.rm=TRUE,names=FALSE),2)
  Q25 <- round(period.apply(ts, spec_period, FUN = quantile,probs=c(0.25),na.rm=TRUE,names=FALSE),2)
  Q50 <- round(period.apply(ts, spec_period, FUN = quantile,probs=c(0.5),na.rm=TRUE,names=FALSE),2)
  Q75 <- round(period.apply(ts, spec_period, FUN = quantile,probs=c(0.75),na.rm=TRUE,names=FALSE),2)
  Q95 <- round(period.apply(ts, spec_period, FUN = quantile,probs=c(0.95),na.rm=TRUE,names=FALSE),2)
  IQR<-abs(Q75-Q25)
  
  period_stats <- cbind.xts(NumofData,NumofMisData,PercOfMissingData,Min,Max,Mean,Var,StDev,Variation,Mom3,Skewness,
                            Kurtosis,LMean,LScale,L3,L4,LVariation,LSkewness,Lkurtosis,Pdr,Q5,Q25,Q50,Q75,Q95,IQR)
  colnames(period_stats) <- c('NumofData','NumofMisData','PercOfMissingData','Min','Max','Mean','Var','StDev',
                             'Variation','Mom3','Skewness','Kurtosis','Lmean','LScale','L3','L4',
                             'LVariation','LSkewness','LKurtosis','Pdr','Q5','Q25','Q50','Q75','Q95','IQR')
  return(period_stats)
}

