#' @title basic_stats
#'
#' @description  This function provides basic statistics and plots of a timeseries in xts format. 
#' It plots the timeseries values, the empirical probability density function (PDF), 
#' the empirical cumulative density function (CDF) and the autocorrelation. Additionally, it will calculate 
#' all major statistics. 
#' \itemize{
#' \item Number of Data - Number of missing data - Percentage of missing Data
#' \item Minimum - Maximum - Average - Variance - Coefficient of Variation - Standard Deviation - Third moment - Skewness 
#' - Kurtosis 
#' \item L-mean, L-scale , 3rd and 4th L-coefficients
#' \item Probability Dry 
#' \item Quantiles 5th , 25th, 50th, 75th, 95th - Interquartile range. 
#' \item Mean value and variance from a dry (defined by zero_threshold argument) to a wet state
#' \item Mean value and variance from a wet to a dry (defined by zero_threshold argument) state
#' \item Mean value and variance from a wet to a wet state
#' \item Transition probability from a wet to a wet state
#' \item Transition probability from a dry to a dry state
#' }
#' With arguments for plotting position of stats table and timeseries plotting period.
#' 
#' @param ts A xts object containing the time series data. 
#' @param pstart Plotting start date. If not included, the first date of the timeseries will be used.
#' Usefull for long timeseries. 
#' @param pend Plotting end date. If not included, the last date of the timeseries will be used.
#' Usefull for long timeseries. 
#' @param show_label A logical value, if TRUE the timeseries title will be plotted. 
#' It is assumed to be the column name. Default is TRUE.
#' @param label_prefix A character value wich contains the prefix of the timeseries title. Default is 'station'
#' @param show_table A logical value, if TRUE a table with the mean value, the standard deviation, the skewness and 
#' the probability dry is plotted together with the timeseries. Default is true. 
#' @param xpos_label The x position of the timeseries label. Takes values from 0 to 1.
#' @param ypos_label The y position of the timeseries label. Takes values from 0 to 1.
#' @param xpos_table The x position of the timeseries table with the basic statistics. Takes values from 0 to 1.
#' @param ypos_table The y position of the timeseries table with the basic statistics. Takes values from 0 to 1.
#' @param nbins The number of bins to split the data from the PDF calculation. Default is 30.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored when computing the statistics. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param title A logical value, if TRUE then a title will be added to the plots. Default is FALSE.
#' 
#' @return A list which contains the combined timeseries, PDF, ECDF and autocorrelation plot and the statistics table. 
#' 
#' @example 
#'file <- "KNMI_Daily.csv"
#'file_path <- file.path(cwd,file)
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path, 
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#' 
#' bstats <- basic_stats(data[,4], pstart = '2002',pend = '2005', show_table = FALSE,
#' show_label = FALSE, ignore_zeros = TRUE)
#' 
#' bstats$plot
#' bstats$stats_table
#' 
#' @export


# Function to plot the basic stats from a single timeseries. 
# Plots the raw ts, PDF, empirical CDF and ACF. Calculates all major statistics. 
# - Number of Data - Number of missing data - Percentage of missing Data - Minimum - Maximum 
# - Average - Variance - Coefficient of Variation - Standard Deviation - Third moment - Skewness 
# - Kurtosis - Probability Dry - Quantiles 5th , 25th, 50th, 75th, 95th - Interquartile range. 
# With arguments for setting zero threshold and ignoring zeros.
# When zeros are ignored ECDF, PDF and stats calculations ignore zeros. 
# With arguments for plotting position of stats table and raw ts plotting period 
# pstart -> Plotting start date 
# pend -> plotting end date
# show_label -> Plot the timeseries title. It is assumed to be the column name 
# label_prefix -> Prefix label e.g. Station, Variable etc.
# show_table -> Plot basic stat table. These include - Average - Coefficient of Variation - Skewness - Probability Dry
# xpos_label, ypos_label etc -> arguments to position the label and table relative the maximum and minimum values of 
# x and y axis. 

basic_stats <- function(ts, pstart = NA, pend = NA, show_label = TRUE, label_prefix = 'Station', show_table = TRUE,
                        xpos_label = 0.1, ypos_label = 0.95,
                         xpos_table = 0.1, ypos_table = 0.15, nbins = 30, ignore_zeros = FALSE, zero_threshold = 0.01){
  if (ignore_zeros == TRUE){
    temp <- ts[ts>zero_threshold,]
  }else{
    temp <- ts
  }
  Mean <- round(mean(temp,na.rm=TRUE),2)
  StDev <- round(sd(temp,na.rm=TRUE),2)
  Variation <- round(StDev/Mean,2)
  Mom3 <- round(sum((temp-Mean)^3,na.rm=T)*(1/(sum(!is.na(temp))-1)),6)
  Skewness<-round(Mom3/StDev^3,2)
  Pdr<-round(mean(ts<=zero_threshold,na.rm=T),2)
  
  stats_basic <- cbind(Metric=c('Mean','Variation','Skewness','Prob. Dry'), Value=c(Mean,Variation,Skewness,Pdr))
  
  if (!is.na(pstart)){start_date <- pstart}else{start_date <- index(ts)[1]}
  if (!is.na(pend)){end_date <- pend}else{end_date <- index(ts)[nrow(ts)]}
  plot_period <- paste(start_date,end_date,sep = '/')
  ts_plot <- ts[plot_period]
  
  data_min <- apply(ts_plot, 2, FUN = min, na.rm = TRUE)
  data_max <- apply(ts_plot, 2, FUN = max, na.rm = TRUE)

  plot_rawts <- autoplot.zoo(ts_plot) + ggtitle('Raw Timeseries')
 
  
  if (show_label == TRUE){
    plot_rawts <- plot_rawts + annotate('label', x=index(ts_plot)[round(xpos_label*nrow(ts_plot),digits = 0)],
                                        y= data_min + (data_max - data_min)*ypos_label, 
                                        label = paste(label_prefix,colnames(ts))) + xlab('Date') 
  }
  
  if (show_table == TRUE){
    plot_rawts <- plot_rawts +  annotation_custom(gridExtra::tableGrob(stats_basic, rows = NULL), xmin = index(ts_plot)[round(xpos_table*nrow(ts_plot),digits = 0)], 
                                                  xmax = index(ts_plot)[round(xpos_table*nrow(ts_plot),digits = 0)], ymin = data_min + (data_max - data_min)*ypos_table, 
                                                  ymax = data_min + (data_max - data_min)*(ypos_table+0.05))
  }
  
  ts_df <- data.frame(x=coredata(ts))
  if (ignore_zeros == TRUE){ ts_df <- as.data.frame(ts_df[ts_df > zero_threshold,])}
  names(ts_df) <- 'X'
  plot_hist <- ggplot(ts_df , aes(x=X)) + 
    geom_histogram(aes(y=..density..), bins = nbins,     # Histogram with density instead of count on y-axis
                   colour='black', fill='white') +
    geom_density(alpha=.2, fill='#FF6666') + # Overlay with transparent density plot 
    labs(x = colnames(ts), y = 'Density') + ggtitle('PDF')
  
  #library(scales)
  
  plot_ecdf <- ggplot(ts_df, aes(X)) + stat_ecdf(geom = 'step') + scale_y_continuous(trans = 'log10', breaks = trans_breaks('log10', function(x) 10^x),
                                                                                     labels = trans_format('log10', math_format(10^.x))) + 
    labs(x = colnames(ts), y = 'P(X<x)') + ggtitle('ECDF')
  
  ts_clean <- ts[!is.na(ts),]
  if (ignore_zeros == TRUE){ts_clean <- ts_clean[which(ts_clean > zero_threshold),]}
  acf_dirty <- acf(ts_clean, lag.max = 10, plot = FALSE)
  acf2df <- data.frame(ACF=acf_dirty$acf, Lag = seq(from = 0, to = 10))
  
  plot_acf <- ggplot(acf2df, aes(x = Lag, y = ACF)) + geom_point() + geom_line() + ggtitle('ACF') #+ ylab(expression(rho_{t,t-1}))
  
  combi_plot <- ggpubr::ggarrange(ggpubr::ggarrange(plot_rawts,plot_hist, widths = c(1.5,1), ncol = 2), 
                                  ggpubr::ggarrange(plot_ecdf,plot_acf, ncol = 2), nrow = 2)
  
  
  NumofData <- length(ts)
  NumofMisData <- length(ts[is.na(ts)])
  PercOfMissingData <- (NumofMisData/NumofData)*100
  Min <- round(min(ts,na.rm=TRUE),2)
  Max <- round(max(ts,na.rm=TRUE),2)
  Mean <- round(mean(ts,na.rm=TRUE),2)
  Var <- round(var(ts,na.rm=TRUE),2)
  StDev <- round(sd(ts,na.rm=TRUE),2)
  Variation <- round(StDev/Mean,2)
  Mom3 <- round(sum((ts-Mean)^3,na.rm=T)*(1/(sum(!is.na(ts))-1)),2)
  Skewness <- round(moments::skewness(ts),2)
  Kurtosis <- round(moments::kurtosis(ts,na.rm=TRUE))
  #L-moments 
  lmom<-lmom::samlmu(coredata(ts), nmom = 4, ratios = FALSE, trim = 0)
  LMean<-round(lmom[1],2)
  LScale<-round(lmom[2],2)
  L3<-round(lmom[3],2)
  L4<-round(lmom[4],2)
  lmom_ratios<-lmom::samlmu(coredata(ts), nmom = 4, ratios = TRUE, trim = 0)
  LVariation<-round(LScale/LMean,2)
  LSkewness<-round(lmom_ratios[3],2)
  Lkurtosis<-round(lmom_ratios[4],2)
  Pdr <- round(mean(ts<=zero_threshold,na.rm=T),2)
  
  # Quantiles of all values
  Q5<-round(quantile(ts,probs=c(0.05),na.rm=TRUE,names=FALSE),2)
  Q25<-round(quantile(ts,probs=c(0.25),na.rm=TRUE,names=FALSE),2)
  Q50<-round(quantile(ts,probs=c(0.50),na.rm=TRUE,names=FALSE),2)  
  Q75<-round(quantile(ts,probs=c(0.75),na.rm=TRUE,names=FALSE),2)
  Q95<-round(quantile(ts,probs=c(0.95),na.rm=TRUE,names=FALSE),2)
  IQR<-abs(Q75-Q25)
  
  # Conditional Statistics and probabilities
  pos<-2:length(ts)
  lagpos<-pos-1
  Ypos<-coredata(ts[pos])
  Ylagpos<-coredata(ts[lagpos])
  
  # Non-zero to Zero
  z<-ifelse(Ypos>zero_threshold & Ylagpos<=zero_threshold,Ypos,NA)
  dem.after.zero<-z[!is.na(z)]
  MeanDAfterZero<-round(mean(dem.after.zero, na.rm=T), 5)
  VarDAfterZero<-round(var(dem.after.zero, na.rm=T), 5)
  # Zero to Non-zero
  z<-ifelse(Ypos<=zero_threshold & Ylagpos>zero_threshold,Ylagpos,NA)
  dem.before.zero<-z[!is.na(z)]
  MeanDBeforeZero<-round(mean(dem.before.zero, na.rm = T),5)
  VarDBeforeZero<-round(var(dem.before.zero, na.rm = T),5)
  # Non-zero to Non-zero
  z<-ifelse(Ypos>zero_threshold & Ylagpos>zero_threshold,Ypos,NA)
  dem.after.dem<-z[!is.na(z)]
  MeanDAfterD<-round(mean(dem.after.dem, na.rm = T),5)
  VarDAfterD<-round(var(dem.after.dem, na.rm = T),5)
  z<-ifelse(Ypos>zero_threshold & Ylagpos>zero_threshold,Ypos,-99)
  ProbDD<-round(length(which(z!=-99 & !is.na(z)))/length(which(!is.na(z))),5)
  # Zero to Zero
  z<-ifelse(Ypos==zero_threshold & Ylagpos==zero_threshold,Ypos,-99)
  ProbNDND <- round(length(which(z!=-99 & !is.na(z)))/length(which(!is.na(z))),5)
  
  stats_table <- data.frame(Value=c(NumofData,NumofMisData,PercOfMissingData,Min,Max,Mean,Var,StDev,Variation,Mom3,Skewness,
                                    Kurtosis,LMean,LScale,L3,L4,LVariation,LSkewness,Lkurtosis,Pdr,Q5,Q25,Q50,Q75,Q95,IQR,
                                    MeanDAfterZero,VarDAfterZero,MeanDBeforeZero,VarDBeforeZero,MeanDAfterD,VarDAfterD,ProbDD,ProbNDND))
  
  rownames(stats_table) <- c('NumofData','NumofMisData','PercOfMissingData','Min','Max','Mean','Var','StDev',
                             'Variation','Mom3','Skewness','Kurtosis','Lmean','LScale','L3','L4',
                            'LVariation','LSkewness','LKurtosis','Pdr','Q5','Q25','Q50','Q75','Q95','IQR',
                            'MeanDAfterZero','VarDAfterZero','MeanDBeforeZero','VarDBeforeZero','MeanDAfterD','VarDAfterD','ProbDD','ProbNDND')
  
  
  if (ignore_zeros == TRUE){
    # Statistics of non-zero values
    MinNonZero<-round(min(ts[ts>zero_threshold],na.rm=TRUE),2)
    MaxNonZero<-round(max(ts[ts>zero_threshold],na.rm=TRUE),2)
    MeanNonZero<-round(mean(ts[ts>zero_threshold],na.rm=T),2)
    VarNonZero<-round(var(ts[ts>zero_threshold],na.rm=T),2)
    StDevNonZero<-round(sd(ts[ts>zero_threshold],na.rm=T),2)
    VariationNonZero <- round(StDevNonZero/MeanNonZero,2)
    Mom3NonZero<-round(sum((ts[ts>zero_threshold]-MeanNonZero)^3,na.rm=T)*(1/(sum(!is.na(ts[ts>zero_threshold]))-1)),2)
    SkewnessNonZero<-round(moments::skewness(ts[ts>zero_threshold],na.rm=TRUE),2)
    KurtosisNonZero <- round(moments::kurtosis(ts[ts>zero_threshold],na.rm=TRUE))
    #L-moments 
    lmom<-lmom::samlmu(coredata(ts[ts>zero_threshold]), nmom = 4, ratios = FALSE, trim = 0)
    LMean<-round(lmom[1],2)
    LScale<-round(lmom[2],2)
    L3<-round(lmom[3],2)
    L4<-round(lmom[4],2)
    lmom_ratios<-lmom::samlmu(coredata(ts[ts>zero_threshold]), nmom = 4, ratios = TRUE, trim = 0)
    LVariation<-round(LScale/LMean,2)
    LSkewness<-round(lmom_ratios[3],2)
    Lkurtosis<-round(lmom_ratios[4],2)
   
    # Quantiles of non-zero values
    Q5NonZero<-round(quantile(ts[ts>zero_threshold],probs=c(0.05),na.rm=TRUE,names=FALSE),5)
    Q25NonZero<-round(quantile(ts[ts>zero_threshold],probs=c(0.25),na.rm=TRUE,names=FALSE),5)
    Q50NonZero<-round(quantile(ts[ts>zero_threshold],probs=c(0.50),na.rm=TRUE,names=FALSE),5)  
    Q75NonZero<-round(quantile(ts[ts>zero_threshold],probs=c(0.75),na.rm=TRUE,names=FALSE),5)
    Q95NonZero<-round(quantile(ts[ts>zero_threshold],probs=c(0.95),na.rm=TRUE,names=FALSE),5)
    IQRNonZero<-abs(Q75NonZero-Q25NonZero)
    
    stats_table <- data.frame(Value=c(NumofData,NumofMisData,PercOfMissingData,MinNonZero,MaxNonZero,MeanNonZero,
                                      VarNonZero,StDevNonZero,VariationNonZero,Mom3NonZero,SkewnessNonZero,
                                      KurtosisNonZero,LMean,LScale,L3,L4,LVariation,LSkewness,Lkurtosis,
                                      Pdr,Q5NonZero,Q25NonZero,Q50NonZero,Q75NonZero,Q95NonZero,IQRNonZero,
                                      MeanDAfterZero,VarDAfterZero,MeanDBeforeZero,VarDBeforeZero,MeanDAfterD,VarDAfterD))
    
    rownames(stats_table) <- c('NumofData','NumofMisData','PercOfMissingData','Min','Max','Mean','Var','StDev',
                               'Variation','Mom3','Skewness','Kurtosis','Lmean','LScale','L3','L4',
                               'LVariation','LSkewness','LKurtosis','Pdr','Q5','Q25','Q50','Q75','Q95','IQR',
                               'MeanDAfterZero','VarDAfterZero','MeanDBeforeZero','VarDBeforeZero','MeanDAfterD','VarDAfterD')
  }
  
  list_out <- list(plot = combi_plot,stats_table = stats_table)
  
  return(list_out)
  
}
