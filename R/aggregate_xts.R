#' @title aggregate_xts
#'
#' @description Aggregates a xts timeseries to different time scales
#' and returns a list of plots and aggregated timeseries.
#'
#' @param ts A xts object containing the time series data.
#' @param periods A vector of periods to aggregate the timeseries to.
#' Accepts 'mins' or 'minutes', 'hours', 'days', 'weeks', 'months', 'quarters', and 'years'.
#' @param FUN The function used to aggregate the data. Defaults to 'mean'.
#' @param period_multiplier Multiplier to aggregate data to a non-standard time scale.
#'  Must be a vector with the same length as periods.
#' @param pstart Plotting start date. If not included, the first date of the timeseries will be used.
#' Usefull for long timeseries.
#' @param pend Plotting end date. If not included, the last date of the timeseries will be used.
#' Usefull for long timeseries.
#' @param season_end Vector of season end dates. Can be used to aggregate into user-defined seasons.
#' Must be in the form of '%m-%d', e.g. '10-15'
#' @param mseason_title Title for the user-defined season.
#'
#' @return A list with the aggregated timeseries for every selected period in xts format and
#' a plot containing the aggregated timeseries.
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#' agg_ts <- aggregate_xts(data[,4], c("months","quarters","years"), FUN = 'sum')
#' smonthly_ts <- agg_ts$list_months$aggregated
#'
#' @export

# Function to aggregate the original timeseries to different scales. Plots the results and saves the aggregated
# timeseries
# pstart, pend <- plot starting and end date for the raw data. Usefull in case of very large ts
# FUN <- Function to aggregate the data. Default is mean. Could be max, sum etc.
# periods <- periods for aggregation. Can be
#“mins” (minutes), “minutes”, “hours”, “days”, “weeks”, “months”, “quarters”, and “years”
# Allows for arbitrary seasons. Need to specify end date for each season.
# season_end <- set of strings in the form of '%m-%d', e.g. '10-15'
# mseason_title <- The title for the manual season.
aggregate_xts <- function(ts, periods = NA, FUN = 'mean', period_multiplier = NA,
                          pstart = NA, pend = NA, season_end = NA, mseason_title = NA){

  if (!is.na(pstart)){start_date <- pstart}else{start_date <- index(ts)[1]}
  if (!is.na(pend)){end_date <- pend}else{end_date <- index(ts)[nrow(ts)]}
  plot_period <- paste(start_date,end_date,sep = '/')

  fraw <- autoplot.zoo(ts[plot_period], geom = 'line') + ggplot2::xlab('Date') +
    ggplot2::ggtitle('Raw Timeseries')

  if (is.na(period_multiplier[1])){
    period_multiplier = rep(1, length(periods))
  }else if(length(period_multiplier) != length(periods)){
    stop('Argument period_multiplier must be a vector with the same length as periods')
  }

  #“secs” (seconds), “seconds”, “mins” (minutes), “minutes”, “hours”, “days”, “weeks”, “months”, “quarters”, and “years”
  suppressWarnings(
    if (!is.na(periods)){
      for (period in periods){
        if (period == 'months'){
          period_title <- 'Monthly'
        }else if(period =='years'){
          period_title <- 'Annual'
        }else if(period == 'quarters'){
          period_title <- 'Quarterly'
        }else if(period == 'weeks'){
          period_title <- 'Weekly'
        }else if(period == 'days'){
          period_title <- 'Daily'
        }else if(period == 'hours'){
          period_title <- 'Hourly'
        }else if(period == 'mins' & period == 'minutes'){
          period_title <- 'Minutes'
        }

        Iper <- which(periods == period)
        spec_period <- endpoints(ts, on = period, k = period_multiplier[Iper])

        period_ts <- period.apply(ts, spec_period, FUN = eval(FUN))

        f <- autoplot.zoo(period_ts, geom = 'line') + ggplot2::xlab('Date') +
          ggplot2::ggtitle(paste(period_multiplier[Iper], period_title,'Timeseries', sep = ' '))
        assign(paste('f', which(periods == period), sep = ''),f)
        assign(paste('list',period, sep = '_'),list(aggregated = period_ts, figure = f))
      }
    }
  )
  suppressWarnings(
    if (!is.na(season_end)){
      for (dt in season_end){
        char_date <- strsplit(dt, '-')
        mmonth <- as.numeric(char_date[[1]][1])
        dday <- as.numeric(char_date[[1]][2])

        Itemp <- which(month(ts)==mmonth & day(ts) == dday)
        if (dt == season_end[1]){
          I <- Itemp
        }else{
          I <- c(I,Itemp)
        }
      }

      spec_period <- sort(I)

      period_ts <- period.apply(ts, spec_period, FUN = eval(FUN))

      f <- autoplot.zoo(period_ts, geom = 'line') + ggplot2::xlab('Date') +
        ggplot2::ggtitle(paste(mseason_title,'Timeseries', sep = ' '))

      if (!is.na(periods)){
        suffix <- length(periods) + 1
      }else{
        suffix <- 1
      }

      assign(paste('f', suffix, sep = ''),f)
      assign(paste('list','manual', sep = '_'),list(aggregated = period_ts, figure = f))

      periods <- c(periods,'manual')
    }
  )

  nperiods <- length(periods)

  if (nperiods==1){
    f_combined <- ggpubr::ggarrange(fraw,f1, nrow = (nperiods + 1), ncol = 1)
    list_out <- list(f_combined, eval(parse(text = paste('list',periods[1],sep = '_'))))
    names(list_out) <- c('Combined_Plot', paste('list',periods[1],sep = '_'))
  }else if(nperiods==2){
    f_combined <- ggpubr::ggarrange(fraw,f1,f2, nrow = (nperiods + 1), ncol = 1)
    list_out <- list(f_combined, eval(parse(text = paste('list',periods[1],sep = '_'))),
                     eval(parse(text = paste('list',periods[2],sep = '_'))))
    names(list_out) <- c('Combined_Plot', paste('list',periods[1],sep = '_'),paste('list',periods[2],sep = '_'))
  }else if(nperiods==3){
    f_combined <- ggpubr::ggarrange(fraw,f1,f2,f3, nrow = (nperiods + 1), ncol = 1)
    list_out <- list(f_combined, eval(parse(text = paste('list',periods[1],sep = '_'))),
                     eval(parse(text = paste('list',periods[2],sep = '_'))),eval(parse(text = paste('list',periods[3],sep = '_'))))
    names(list_out) <- c('Combined_Plot', paste('list',periods[1],sep = '_'),paste('list',periods[2],sep = '_'),
                         paste('list',periods[3],sep = '_'))
  }else if(nperiods==4){
    f_combined <- ggpubr::ggarrange(fraw,f1,f2,f3,f4, nrow = (nperiods + 1), ncol = 1)
    list_out <- list(f_combined, eval(parse(text = paste('list',periods[1],sep = '_'))),
                     eval(parse(text = paste('list',periods[2],sep = '_'))),eval(parse(text = paste('list',periods[3],sep = '_'))),
                     eval(parse(text = paste('list',periods[4],sep = '_'))))
    names(list_out) <- c('Combined_Plot', paste('list',periods[1],sep = '_'),paste('list',periods[2],sep = '_'),
                         paste('list',periods[3],sep = '_'),paste('list',periods[4],sep = '_'))
  }else if(nperiods==5){
    f_combined <- ggpubr::ggarrange(fraw,f1,f2,f3,f4,f5, nrow = (nperiods + 1), ncol = 1)
    list_out <- list(f_combined, eval(parse(text = paste('list',periods[1],sep = '_'))),
                     eval(parse(text = paste('list',periods[2],sep = '_'))),eval(parse(text = paste('list',periods[3],sep = '_'))),
                     eval(parse(text = paste('list',periods[4],sep = '_'))), eval(parse(text = paste('list',periods[5],sep = '_'))))
    names(list_out) <- c('Combined_Plot', paste('list',periods[1],sep = '_'),paste('list',periods[2],sep = '_'),
                         paste('list',periods[3],sep = '_'),paste('list',periods[4],sep = '_'),paste('list',periods[5],sep = '_'))
  }else if(nperiods==6){
    f_combined <- ggpubr::ggarrange(fraw,f1,f2,f3,f4,f5,f6, nrow = (nperiods + 1), ncol = 1)
    list_out <- list(f_combined, eval(parse(text = paste('list',periods[1],sep = '_'))),
                     eval(parse(text = paste('list',periods[2],sep = '_'))),eval(parse(text = paste('list',periods[3],sep = '_'))),
                     eval(parse(text = paste('list',periods[4],sep = '_'))), eval(parse(text = paste('list',periods[5],sep = '_'))),
                     eval(parse(text = paste('list',periods[6],sep = '_'))))
    names(list_out) <- c('Combined_Plot', paste('list',periods[1],sep = '_'),paste('list',periods[2],sep = '_'),
                         paste('list',periods[3],sep = '_'),paste('list',periods[4],sep = '_'),paste('list',periods[5],sep = '_'),
                         paste('list',periods[6],sep = '_'))
  }

  return(list_out)
}
