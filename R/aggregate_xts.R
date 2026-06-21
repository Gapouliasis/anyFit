#' Aggregate an xts time series to coarser temporal scales
#'
#' @description
#' Aggregates an xts object to one or more coarser temporal scales using
#' \code{\link[xts]{endpoints}} and \code{\link[xts]{period.apply}}. Supports
#' user-defined seasonal aggregation via \code{season_end} dates and non-integer
#' period multipliers. Returns both the aggregated series and a combined
#' \code{patchwork} plot.
#'
#' @param ts An xts object containing the time series data.
#' @param periods Character vector of aggregation periods. Accepts \code{"mins"},
#'   \code{"minutes"}, \code{"hours"}, \code{"days"}, \code{"weeks"},
#'   \code{"months"}, \code{"quarters"}, and \code{"years"}.
#' @param FUN Function to aggregate the data. Defaults to \code{"mean"}.
#' @param period_multiplier Numeric vector of multipliers for non-standard
#'   aggregation intervals. Must have the same length as \code{periods}.
#'   Defaults to 1 for each period.
#' @param pstart Plot start date (character or date). If \code{NA}, the first
#'   date of the series is used.
#' @param pend Plot end date (character or date). If \code{NA}, the last date
#'   of the series is used.
#' @param season_end Character vector of season end dates in \code{"\%m-\%d"}
#'   format, e.g. \code{"10-15"}.
#' @param mseason_title Title for the user-defined season plot.
#' @param ... Additional arguments passed to \code{FUN}.
#'
#' @return A named list. The first element \code{Combined_Plot} is a
#'   \code{patchwork} object combining all period plots. Subsequent elements
#'   are named \code{list_<period>}, each a sub-list containing the
#'   \code{aggregated} xts and its \code{figure}.
#'
#' @examples
#' # Synthetic daily precipitation series
#' set.seed(123)
#' n <- 1000
#' x <- xts::xts(rgamma(n, shape = 0.6, scale = 5),
#'               order.by = seq(as.Date("2000-01-01"), by = "day", length.out = n))
#'
#' # Aggregate to monthly sums
#' agg <- aggregate_xts(x, c("months", "years"), FUN = "sum")
#' agg$Combined_Plot
#'
#' # Custom seasonal aggregation
#' agg_seas <- aggregate_xts(x, periods = NA, FUN = "sum",
#'                            season_end = "09-30",
#'                            mseason_title = "Custom season")
#'
#' @importFrom lubridate day
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
                          pstart = NA, pend = NA, season_end = NA, mseason_title = NA, ...){

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

  # if (is.na(period_multiplier)){
  #   period_multiplier = rep(1,length(periods))
  # }

  plot_list = list(fraw)
  agg_list = list()
  #“secs” (seconds), “seconds”, “mins” (minutes), “minutes”, “hours”, “days”, “weeks”, “months”, “quarters”, and “years”
  suppressWarnings(
    if (!is.na(periods)[1]){
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

        period_ts <- period.apply(ts, spec_period, FUN = eval(FUN))#,...)

        f <- autoplot.zoo(period_ts, geom = 'line') + ggplot2::xlab('Date') +
          ggplot2::ggtitle(paste(period_multiplier[Iper], period_title,'Timeseries', sep = ' '))
        plot_list = c(plot_list, list(f))

        agg_list = c(agg_list, list(list(aggregated = period_ts, figure = f)))
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

      plot_list = c(plot_list, list(f))

      agg_list = c(agg_list, list(aggregated = period_ts, figure = f))

      periods <- c(periods,'manual')
    }
  )

  nperiods <- length(periods)

  f_combined = patchwork::wrap_plots(plotlist = plot_list, nrow = (nperiods + 1), ncol = 1)
  list_out = c(list(f_combined), agg_list)
  names(list_out) = c('Combined_Plot', paste('list',periods,sep = '_'))

  return(list_out)
}
