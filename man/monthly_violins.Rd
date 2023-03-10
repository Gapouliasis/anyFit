% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monthly_violins.R
\name{monthly_violins}
\alias{monthly_violins}
\title{monthly_violins}
\usage{
monthly_violins(
  ts,
  palette = "Set3",
  ignore_zeros = FALSE,
  zero_threshold = 0.01
)
}
\arguments{
\item{ts}{A xts object containing the time series data}

\item{palette}{The color palette to use. Default is 'Set3'.}

\item{ignore_zeros}{A logical value, if TRUE zeros will be ignored. Default is FALSE.}

\item{zero_threshold}{The threshold below which values are considered zero. Default is 0.01.}
}
\value{
A ggplot2 object.
}
\description{
Plots a violin plot of monthly values
}
\examples{
file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
time_zone <- "UTC"
time_step <- "1 day"

data <- delim2xts(file_path = file_path,
                 time_zone = "UTC", delim = " ", time_step = time_step)

violins <- monthly_violins(data[,4], palette = 'Set3', ignore_zeros = TRUE)

}
