#' @title monthly_stats
#'
#' @description  This function provides monthly statistics in aggregated and original (base) timescale 
#' and plots of a time series. Statistics ploted are average, standard deviation, skewness. For the base scale
#' the probability dry is calculated. For the monthly scale the month to month correlations are calculated. 
#' Month to month correlations are calculated backwards. i.e January == cor(December, January). Statistics calculated 
#' and provided on a table are the number of data points, 
#' number of missing data, percentage of missing data, min, max, mean, variance, 
#' standard deviation, variation, 3rd moment, skewness, kurtosis and probability dry. 
#'
#' @param ts A xts object containing the time series data. 
#' @param aggregated A logical value, if TRUE the statistics are also given on the aggregated monthly scale
#'  of the time series. Default is FALSE.
#' @param FUN The aggregation function applied to the time series. Default is 'mean'.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored when computing the statistics. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param title A logical value, if TRUE then a title will be added to the plots. Default is FALSE.
#' @return If aggregated is TRUE, a list containing the aggregated statistics and the plots is returned. 
#' If aggregated is FALSE, a list containing the base scale statistics and the plots is returned.
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
#' monthly_sts <- monthly_stats(data[,4],base_scale = TRUE, FUN = "sum")
#' 
#' monthly_sts$fbase
#' monthly_sts$base_list$July$raw_stats
#' 
#' @export

#Function to calculate the seasonal statistics in aggregated and original (base) timescale. 
#Correlations are calculated backwards. i.e January == cor(December, January)
monthly_stats <- function(ts,aggregated = FALSE, FUN = 'mean', ignore_zeros = FALSE, zero_threshold = 0.01, title = FALSE){
# base_scale <- TRUE
# FUN = 'mean'
#ts <- ts[which(!is.na(ts)),]
  monthly_ts <- apply.monthly(ts, FUN = eval(FUN))
  agg_means <- data.frame(Month = 1:12, Value = rep(0,12))
  agg_stdev <- data.frame(Month = 1:12, Value = rep(0,12))
  agg_skews <- data.frame(Month = 1:12, Value = rep(0,12))
  agg_correls <- data.frame(Month = 1:12, Value = rep(0,12))
  for (i in 1:12){
    I <- which(month(monthly_ts) == i)
    stats <- basic_stats(monthly_ts[I],ignore_zeros = ignore_zeros, zero_threshold = zero_threshold)
    temp_stats <- stats$stats_table
    colnames(temp_stats) <- month.name[i]
    if (i == 1){
      agg_stats <- temp_stats
    }else{
      agg_stats <- cbind(agg_stats, temp_stats)
    }
    agg_means[i,2] <- as.numeric(stats$stats_table[6,])
    agg_stdev[i,2] <- as.numeric(stats$stats_table[8,])
    agg_skews[i,2] <- as.numeric(stats$stats_table[11,])
    lag_date <- as.POSIXct(format(index(monthly_ts[I-1]) %m+% months(1), format = '%m/%Y'), format = '%d/%Y', tz = time_zone)
    lag_xts <- xts(zoo::coredata(monthly_ts[I-1]), order.by = lag_date)
    t0_date <- as.POSIXct(format(index(monthly_ts[I]), format = '%m/%Y'), format = '%d/%Y', tz = time_zone)
    t0_xts <-  xts(zoo::coredata(monthly_ts[I]), order.by = t0_date)
    alligned_xts <- merge.xts(t0_xts, lag_xts,join = 'inner')
    agg_correls[i,2] <- round(stats::cor(zoo::coredata(alligned_xts[,1]),zoo::coredata(alligned_xts[,2]), method = 'pearson', use = 'complete.obs'), digits = 2)
    
    # if (i==1){
    #   agg_list <- list(list(agg_stats = agg_stats, agg_quantiles = agg_quantiles)) 
    # }else{
    #   agg_list <- c(agg_list,list(list(agg_stats = agg_stats, agg_quantiles = agg_quantiles)))
    # }
  }
  
  #names(agg_list) <- month.name[seq(1,12)]
  
  agg_means$Month <- factor(x = month.name[agg_means$Month], levels = month.name)
  plt_mean <- ggplot() + geom_bar(data = agg_means, aes(x = Month, y= Value, fill = as.factor(Value)), color='black', stat = 'identity') + 
    scale_fill_grey(start=0.9,end =0) + theme(legend.position = 'none') + ylab('Average') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#+ ggtitle('Monthly Average')
  
  agg_stdev$Month <- factor(x = month.name[agg_means$Month], levels = month.name)
  plt_stdev <- ggplot() + geom_bar(data = agg_stdev, aes(x = Month, y= Value, fill = as.factor(Value)), color='black', stat = 'identity') + 
    scale_fill_grey(start=0.9,end =0) + theme(legend.position = 'none') + ylab('St. Devation') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#+ ggtitle('Monthly Standard Devation')
  
  agg_skews$Month <- factor(x = month.name[agg_means$Month], levels = month.name)
  plt_skews <- ggplot() + geom_bar(data = agg_skews, aes(x = Month, y= Value, fill = as.factor(Value)), color='black', stat = 'identity') + 
    scale_fill_grey(start=0.9,end =0) + theme(legend.position = 'none') + ylab('Skewness') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#+ ggtitle('Monthly Skewness')
  
  agg_correls$Month <- factor(x = month.name[agg_means$Month], levels = month.name)
  plt_correls <- ggplot() + geom_bar(data = agg_correls, aes(x = Month, y= Value, fill = as.factor(Value)), color='black', stat = 'identity') + 
    scale_fill_grey(start=0.9,end =0) + theme(legend.position = 'none') + ylab('Correlation') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#+ ggtitle('Monthly Correlation')
  
  faggre <- ggpubr::ggarrange(plt_mean, plt_stdev, plt_skews, plt_correls, nrow = 2, ncol = 2)
  if (title == TRUE){
    faggre <- ggpubr::annotate_figure(faggre, top = ggpubr::text_grob('Monthly Scale Statistics', 
                                                                      face = 'bold', size = 14),fig.lab.pos = c('top.left'))
  }

  
  if (aggregated == TRUE){
    raw_means <- data.frame(Month = 1:12, Value = rep(0,12))
    raw_stdev <- data.frame(Month = 1:12, Value = rep(0,12))
    raw_skews <- data.frame(Month = 1:12, Value = rep(0,12))
    raw_correls <- data.frame(Month = 1:12, Value = rep(0,12))
    raw_pdr <- data.frame(Month = 1:12, Value = rep(0,12))
    for (i in 1:12){
      I <- which(month(ts) == i)
      stats <- basic_stats(ts[I],ignore_zeros = ignore_zeros, zero_threshold = zero_threshold)
      temp_stats <- stats$stats_table
      colnames(temp_stats) <- month.name[i]
      if (i == 1){
        raw_stats <- temp_stats
      }else{
        raw_stats <- cbind(raw_stats, temp_stats)
      }
      raw_means[i,2] <- as.numeric(stats$stats_table[6,])
      raw_stdev[i,2] <- as.numeric(stats$stats_table[8,])
      raw_skews[i,2] <- as.numeric(stats$stats_table[11,])
      raw_pdr[i,2] <- as.numeric(stats$stats_table[13,])
      lag_date <- as.POSIXct(format(index(ts[I-1]) %m+% months(1), format = '%m/%Y'), format = '%d/%Y', tz = time_zone)
      lag_xts <- xts(zoo::coredata(ts[I-1]), order.by = lag_date)
      t0_date <- as.POSIXct(format(index(ts[I]), format = '%m/%Y'), format = '%d/%Y', tz = time_zone)
      t0_xts <-  xts(zoo::coredata(ts[I]), order.by = t0_date)
      alligned_xts <- merge.xts(t0_xts, lag_xts,join = 'inner')
      raw_correls[i,2] <- round(stats::cor(zoo::coredata(alligned_xts[,1]),zoo::coredata(alligned_xts[,2]), method = 'pearson', use = 'complete.obs'), digits = 2)
      
      # if (i==1){
      #   raw_list <- list(list(raw_stats = raw_stats, raw_quantiles = raw_quantiles))
      # }else{
      #   raw_list <- c(raw_list,list(list(raw_stats = raw_stats, raw_quantiles = raw_quantiles)))
      # }
    }
    
    #names(raw_list) <- month.name[seq(1,12)]
    
    raw_means$Month <- factor(x = month.name[agg_means$Month], levels = month.name)
    pltr_mean <- ggplot() + geom_bar(data = raw_means, aes(x = Month, y= Value, fill = as.factor(Value)), color='black', stat = 'identity') + 
      scale_fill_grey(start=0.9,end =0) + theme(legend.position = 'none') + ylab('Average') +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#+ ggtitle('Monthly Average - Base Scale')
    
    raw_stdev$Month <- factor(x = month.name[agg_means$Month], levels = month.name)
    pltr_stdev <- ggplot() + geom_bar(data = raw_stdev, aes(x = Month, y= Value, fill = as.factor(Value)), color='black', stat = 'identity') + 
      scale_fill_grey(start=0.9,end =0) + theme(legend.position = 'none') + ylab('St. Devation') +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#+ ggtitle('Monthly Standard Devation - Base Scale')
    
    raw_skews$Month <- factor(x = month.name[agg_means$Month], levels = month.name)
    pltr_skews <- ggplot() + geom_bar(data = raw_skews, aes(x = Month, y= Value, fill = as.factor(Value)), color='black', stat = 'identity') + 
      scale_fill_grey(start=0.9,end =0) + theme(legend.position = 'none') + ylab('Skewness') +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#+ ggtitle('Monthly Skewness - Base Scale')
    
    raw_correls$Month <- factor(x = month.name[agg_means$Month], levels = month.name)
    pltr_correls <- ggplot() + geom_bar(data = raw_correls, aes(x = Month, y= Value, fill = as.factor(Value)), color='black', stat = 'identity') + 
      scale_fill_grey(start=0.9,end =0) + theme(legend.position = 'none') + ylab('Correlation') +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#+ ggtitle('Monthly Correlation - Base Scale')
    
    raw_pdr$Month <- factor(x = month.name[agg_means$Month], levels = month.name)
    pltr_pdr <- ggplot() + geom_bar(data = raw_pdr, aes(x = Month, y= Value, fill = as.factor(Value)), color='black', stat = 'identity') + 
      scale_fill_grey(start=0.9,end =0) + theme(legend.position = 'none') + ylab('Probability Dry') +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#+ ggtitle('Monthly Probability Dry - Base Scale')
    
    fraw <- ggpubr::ggarrange(pltr_mean, pltr_stdev, pltr_skews, pltr_pdr, nrow = 2, ncol = 2)
    if (title == TRUE){
      fraw <- ggpubr::annotate_figure(fraw, top = ggpubr::text_grob('Monthly Statistics on original scale', 
                                                                    face = 'bold', size = 14),fig.lab.pos = c('top.left'))
    }

  }
  
  if (aggregated == TRUE){
    list_out <- list(agg_stats = agg_stats, faggre = faggre, lag1 = agg_correls)
  }else{
    list_out <- list(base_stats = raw_stats, fbase = fraw)
  }
  
  return(list_out)
}




