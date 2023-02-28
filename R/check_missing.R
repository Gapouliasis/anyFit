#' @title check_missing
#'
#' @description This function takes in time series data and one or more periods and checks for missing values.
#' It returns a data frame with the percentage of missing values for each period and a plot if set to TRUE.
#' The plot shows the percentage of missing values for each period.
#' If the period is set to "months" and "group_months" is set to TRUE,
#' the function returns also a plot of the percentage of missing values grouped by month.
#'
#' @param data A xts object or matrix containing the timeseries data
#' @param periods A character vector containing one or more of the following:
#' "mins" or "minutes", "hours", "days", "weeks", "months", "quarters", and "years".
#' @param plot A logical value indicating whether to plot the missing values or not
#' @param group_months A logical value indicating whether to group the monthly missing values or not
#'
#' @return A list containing a data frame of the percentage of missing values for each period, and optionally a plot if set to TRUE
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#' missing <- check_missing(data[,seq(1,4)], c("months","years"))
#'
#' missing$list_years$figure
#'
#'
#' @export
#'

#function that takes an xts argument (i.e. from delim2xts) and finds the overall percentage of missing values
#and the missing values per period prescribed (i.e. month, year etc).
#Returns the timeseries with the percentage of missing values for the selected period(s) and its plot
check_missing <-function(data, periods, plot = TRUE ,group_months = FALSE){
  NA_table <- apply(data,2, FUN = is.na)
  prct_missing <- as.data.frame(apply(NA_table,2, FUN = sum)/nrow(NA_table) * 100)
  colnames(prct_missing) <- 'Percentage Missing'

  NA_xts <- xts(NA_table, order.by = index(data))

  column_sum <- function(x){
    apply(x, 2, FUN = sum)
  }

  #“secs” (seconds), “seconds”, “mins” (minutes), “minutes”, “hours”, “days”, “weeks”, “months”, “quarters”, and “years”
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


    f <- autoplot(period_prct, geom = 'bar') + ggplot2::xlab('Date') +
      ggplot2::ggtitle(paste('Pecentage of missing values per',period_title))

    if (period == 'months' & group_months == TRUE){
      monthly_missing <- matrix(data = NA,12, ncol(NA_xts))
      colnames(monthly_missing) <- colnames(NA_xts)
      for (i in 1:12){
        monthly_ts <- NA_xts[month(index(NA_xts)) == i,]
        monthly_missing[i,] <- apply(monthly_ts, MARGIN = 2, FUN = sum)/nrow(monthly_ts) * 100
      }
      rownames(monthly_missing) <- month.name
      period_prct = list(monthly_prct = period_prct,grouped_months = monthly_missing)
    }

    if (plot == TRUE){
      assign(paste('list',period, sep = '_'),list(prct_missing = period_prct, figure = f))
    }else{
      assign(paste('list',period, sep = '_'),list(prct_missing = period_prct))
    }

  }

  if (length(periods)==1){
    list_out <- list(prct_missing, eval(parse(text = paste('list',periods[1],sep = '_'))))
    names(list_out) <- c('prct_missing', paste('list',periods[1],sep = '_'))
  }else if(length(periods)==2){
    list_out <- list(prct_missing, eval(parse(text = paste('list',periods[1],sep = '_'))),
                     eval(parse(text = paste('list',periods[2],sep = '_'))))
    names(list_out) <- c('prct_missing', paste('list',periods[1],sep = '_'),paste('list',periods[2],sep = '_'))
  }else if(length(periods)==3){
    list_out <- list(prct_missing, eval(parse(text = paste('list',periods[1],sep = '_'))),
                     eval(parse(text = paste('list',periods[2],sep = '_'))),eval(parse(text = paste('list',periods[3],sep = '_'))))
    names(list_out) <- c('prct_missing', paste('list',periods[1],sep = '_'),paste('list',periods[2],sep = '_'),
                         paste('list',periods[3],sep = '_'))
  }else if(length(periods)==4){
    list_out <- list(prct_missing, eval(parse(text = paste('list',periods[1],sep = '_'))),
                     eval(parse(text = paste('list',periods[2],sep = '_'))),eval(parse(text = paste('list',periods[3],sep = '_'))),
                     eval(parse(text = paste('list',periods[4],sep = '_'))))
    names(list_out) <- c('prct_missing', paste('list',periods[1],sep = '_'),paste('list',periods[2],sep = '_'),
                         paste('list',periods[3],sep = '_'),paste('list',periods[4],sep = '_'))
  }else if(length(periods)==5){
    list_out <- list(prct_missing, eval(parse(text = paste('list',periods[1],sep = '_'))),
                     eval(parse(text = paste('list',periods[2],sep = '_'))),eval(parse(text = paste('list',periods[3],sep = '_'))),
                     eval(parse(text = paste('list',periods[4],sep = '_'))), eval(parse(text = paste('list',periods[5],sep = '_'))))
    names(list_out) <- c('prct_missing', paste('list',periods[1],sep = '_'),paste('list',periods[2],sep = '_'),
                         paste('list',periods[3],sep = '_'),paste('list',periods[4],sep = '_'),paste('list',periods[5],sep = '_'))
  }else if(length(periods)==6){
    list_out <- list(prct_missing, eval(parse(text = paste('list',periods[1],sep = '_'))),
                     eval(parse(text = paste('list',periods[2],sep = '_'))),eval(parse(text = paste('list',periods[3],sep = '_'))),
                     eval(parse(text = paste('list',periods[4],sep = '_'))), eval(parse(text = paste('list',periods[5],sep = '_'))),
                     eval(parse(text = paste('list',periods[6],sep = '_'))))
    names(list_out) <- c('prct_missing', paste('list',periods[1],sep = '_'),paste('list',periods[2],sep = '_'),
                         paste('list',periods[3],sep = '_'),paste('list',periods[4],sep = '_'),paste('list',periods[5],sep = '_'),
                         paste('list',periods[6],sep = '_'))
  }
  return(list_out)
}
