% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitlm_multi.R
\name{fitlm_multi}
\alias{fitlm_multi}
\title{fitlm_multi}
\usage{
fitlm_multi(ts, candidates, ignore_zeros = FALSE, zero_threshold = 0.01)
}
\arguments{
\item{ts}{A xts object containing the time series data.}

\item{candidates}{A list of distribution to fit.}

\item{ignore_zeros}{A logical value, if TRUE zeros will be ignored. Default is FALSE.}

\item{zero_threshold}{The threshold below which values are considered zero. Default is 0.01.}
}
\value{
A list containing the fitted parameters, Goodness-of-Fit Summary, and diagnostic plots.
}
\description{
This function fits a list of candidate distributions using the L-moments method
to a timeseries in xts format. Additionally to the list of the fitted parameters,
goodness-of-fit metric, PP and QQ plots.
}
\examples{
file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
time_zone <- "UTC"
time_step <- "1 day"

data <- delim2xts(file_path = file_path,
                 time_zone = "UTC", delim = " ", time_step = time_step)

candidates <- list('exp','expweibull', 'gamma3')
fits <- fitlm_multi(data['2010',4],candidates = candidates, ignore_zeros = TRUE)

fits$diagnostics

}
