#' Check for missing values in time series data by period
#'
#' @description
#' Computes the percentage of missing (\code{NA}) values in an xts object,
#' both overall and aggregated by one or more time periods. Assumes the timeseries have 
#' been loaded using \code{delim2xts} function so that missing dates are filled with NA. 
#' If this is not the case, or the timestep is not strict then the resulting missing values
#' will be wrong. Optionally groups
#' monthly missingness counts by calendar month for seasonal diagnostics.
#' Returns a data frame of overall missing percentages and, for each period,
#' the period-level percentages and an optional bar plot.
#'
#' @param data An xts object or matrix containing the time series data.
#' @param periods Character vector of aggregation periods. Accepts \code{"mins"},
#'   \code{"minutes"}, \code{"hours"}, \code{"days"}, \code{"weeks"},
#'   \code{"months"}, \code{"quarters"}, and \code{"years"}.
#' @param plot Logical. If \code{TRUE}, bar plots of missing percentages per
#'   period are included in the output. Default is \code{TRUE}.
#' @param group_months Logical. If \code{TRUE} and \code{"months"} is included
#'   in \code{periods}, the function additionally returns a matrix of missing
#'   percentages grouped by calendar month (rows) and series column (columns).
#'   Default is \code{FALSE}.
#'
#' @return A named list. The first element \code{prct_missing} is a data frame
#'   of overall missing percentages (one row per series column). Subsequent
#'   elements are named \code{list_<period>}, each a sub-list containing
#'   \code{prct_missing} (the period-aggregated missing percentages) and, if
#'   \code{plot = TRUE}, \code{figure} (a ggplot bar plot).
#'
#' @examples
#' # Synthetic daily data with injected NAs
#' set.seed(123)
#' n <- 1000
#' x <- xts::xts(cbind(
#'   S1 = rgamma(n, shape = 0.6, scale = 5),
#'   S2 = rgamma(n, shape = 0.7, scale = 4)),
#'   order.by = seq(as.Date("2000-01-01"), by = "day", length.out = n))
#' # Inject ~5% missing values
#' x[sample(n, 50), 1] <- NA
#' x[sample(n, 50), 2] <- NA
#'
#' # Check missingness by month and year
#' miss <- check_missing(x, c("months", "years"))
#' miss$prct_missing
#' miss$list_years$figure
#'
#' # With monthly grouping
#' miss <- check_missing(x, c("months"), group_months = TRUE)
#' miss$list_months$prct_missing$grouped_months
#'
#' @importFrom zoo autoplot.zoo
#' @importFrom xts xts endpoints period.apply
#' @importFrom zoo index coredata
#' @importFrom matrixStats colSums2
#'
#' @export
#'

#function that takes an xts argument (i.e. from delim2xts) and finds the overall percentage of missing values
#and the missing values per period prescribed (i.e. month, year etc).
#Returns the timeseries with the percentage of missing values for the selected period(s) and its plot
check_missing <-function(data, periods, plot = TRUE ,group_months = FALSE){
  NA_table <- is.na(zoo::coredata(data))
  na_counts <- matrixStats::colSums2(NA_table)
  names(na_counts) <- colnames(NA_table)
  prct_missing <- as.data.frame(na_counts/nrow(NA_table) * 100)
  colnames(prct_missing) <- 'Percentage Missing'

  NA_xts <- xts(NA_table, order.by = index(data))

  column_sum <- function(x){
    matrixStats::colSums2(x)
  }

  temp_list = list()

  #"secs" (seconds), "seconds", "mins" (minutes), "minutes", "hours", "days", "weeks", "months", "quarters", and "years"
  for (period in periods){
    if (period == 'months'){
      period_title <- 'Month'
    }else if(period =='years'){
      period_title <- 'Years'
    }else if(period == 'quarters'){
      period_title <- 'Quarter'
    }else if(period == 'weeks'){
      period_title <- 'Week'
    }else if(period == 'days'){
      period_title <- 'Day'
    }else if(period == 'mins' & period == 'minutes'){
      period_title <- 'Minutes'
    }

    spec_period <- endpoints(NA_xts, on = period)

    nperiod <- diff(spec_period, lag = 1)

    period_prct <- period.apply(NA_xts, spec_period, FUN = column_sum)/nperiod * 100
    #prct_df <- fortify.zoo(period_prct) ggplot2::ggplot(data = temp, aes(x=Time,y=Value)) + geom_bar(stat='identity')
    #temp <- prct_df[,c(1,2)]
    #colnames(temp) <- c('Time','Value')


    f <- autoplot.zoo(period_prct, ts.geom = 'bar') + ggplot2::xlab('Date') +
      ggplot2::ggtitle(paste('Pecentage of missing values per',period_title))

    if (period == 'months' & group_months == TRUE){
      monthly_missing <- matrix(data = NA,12, ncol(NA_xts))
      colnames(monthly_missing) <- colnames(NA_xts)
      for (i in 1:12){
        monthly_ts <- NA_xts[month(index(NA_xts)) == i,]
        monthly_missing[i,] <- matrixStats::colSums2(monthly_ts)/nrow(monthly_ts) * 100
      }
      rownames(monthly_missing) <- month.name
      period_prct = list(monthly_prct = period_prct,grouped_months = monthly_missing)
    }

    if (plot == TRUE){
      temp_list = c(temp_list, list(list(prct_missing = period_prct, figure = f)))
    }else{
      temp_list = c(temp_list, list(list(prct_missing = period_prct)))
    }

  }

  list_out <- c(prct_missing, temp_list)
  names(list_out) <- c('prct_missing', paste('list',periods,sep = '_'))

  return(list_out)
}
