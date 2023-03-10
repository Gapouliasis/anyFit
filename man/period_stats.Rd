% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/period_stats.R
\name{period_stats}
\alias{period_stats}
\title{period_stats}
\usage{
period_stats(ts, period = NA, period_multiplier = 1)
}
\arguments{
\item{ts}{a xts object containing the time series data}

\item{period}{a period for which to calculate statistics. Default is NA}

\item{period_multiplier}{a multiplier to define arbitrary periods based on the period argument.
E.g. 2 months, 3 years etc. Default is 1.}
}
\value{
a xts object containing the period-based statistics for the time series data.
}
\description{
Function to calculate period based basic statistics. Specifically, the number of data points,
number of missing data, percentage of missing data, min, max, mean, variance,
standard deviation, variation, 3rd moment, skewness, kurtosis, l-statistics, i.e. mean, scale,
3rd and 4th order l-moments, l-moment variation, l-moment skewness, l-moment kurtosis,
quantiles of 5, 25, 50, 75 and 95, and the inter quartile range.
}
\examples{
file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
time_zone <- "UTC"
time_step <- "1 day"

data <- delim2xts(file_path = file_path,
                 time_zone = "UTC", delim = " ", time_step = time_step)

period_stats(ts, period = "months", period_multiplier = 3)

}
