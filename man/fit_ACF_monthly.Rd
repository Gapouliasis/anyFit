% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_ACF_monthly.R
\name{fit_ACF_monthly}
\alias{fit_ACF_monthly}
\title{fit_ACF_monthly}
\usage{
fit_ACF_monthly(ts, lag, type = list("CAS", "HK", "SRD"), nrow = 4, ncol = 3)
}
\arguments{
\item{ts}{A xts object containing the time series data.}

\item{type}{A list of character strings, containing the type of ACF to fit.
Options include: \code{'CAS'}, \code{'HK'}, and \code{'SRD'}.}

\item{nrow}{Number of rows for plotting.}

\item{ncol}{Number of columns for plotting.}

\item{lag_max}{Maximum lag to use in fitting.}

\item{ignore_zeros}{A logical value, if TRUE zeros will be ignored. Default is FALSE.}

\item{zero_threshold}{The threshold below which values are considered zero. Default is 0.01.}
}
\value{
A list containing a vector of monthly fitted parameters, a data frame of
fitted ACF values, and a plot of the fitted ACFs.
}
\description{
Function to fit theoretical Autocorrelation Functions (ACFs) to monthly values of timeseries data.
Three functions are available, the Cauchy-type autocorrelation structure,
the Hurst - Kolmogorov and the Short Range dependence structure.
}
\examples{
file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
time_zone <- "UTC"
time_step <- "1 day"

data <- delim2xts(file_path = file_path,
                 time_zone = "UTC", delim = " ", time_step = time_step)

monthly_acfs <- fit_ACF_monthly(data[,4], lag = 10)
monthly_acfs$ACF_monthly_plot

}
