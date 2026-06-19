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
#' @param plot A logical value, if TRUE the timeseries/PDF/ECDF/ACF panel is built. Default is FALSE.
#' For a multi-column input one panel is built per column.
#'
#' @return A list with two elements: \code{stats_table}, a data.frame of statistics
#' (rows are the statistics, columns are the input series - named after the input columns
#' when more than one is supplied, otherwise \code{"Value"}); and \code{plot}, which is
#' \code{NULL} when \code{plot = FALSE}, a single combined plot for a single-column input,
#' or a named list of one combined plot per column for a multi-column input.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density labs ggtitle stat_ecdf scale_y_continuous geom_point geom_line annotate annotation_custom xlab theme element_text after_stat
#' @importFrom patchwork plot_layout
#' @importFrom scales trans_breaks trans_format math_format
#' @importFrom zoo autoplot.zoo coredata index
#' @importFrom matrixStats colMins colMaxs colMeans2 colVars colSds colSums2 colQuantiles
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#' bstats <- basic_stats(data[,4], pstart = '2002',pend = '2005', show_table = FALSE,
#' show_label = FALSE, ignore_zeros = TRUE, plot = TRUE)
#'
#' bstats$plot
#' bstats$stats_table
#'
#' # Wide input: one column of statistics per series, plus one plot per series
#' wstats <- basic_stats(data[, 1:4], plot = TRUE)
#' wstats$stats_table          # columns named after the input series
#' wstats$plot[[1]]            # plot for the first series
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
                         xpos_table = 0.1, ypos_table = 0.15, nbins = 30, ignore_zeros = FALSE, zero_threshold = 0.01,
                        plot = FALSE){

  m <- zoo::coredata(ts)
  if (is.null(dim(m))) m <- matrix(m, ncol = 1)
  k <- ncol(m)
  vnames <- colnames(ts)
  if (is.null(vnames)) vnames <- if (k == 1) 'Value' else paste0('V', seq_len(k))

  core <- .basic_stats_core(m, zero_threshold = zero_threshold, ignore_zeros = ignore_zeros)
  stats_table <- as.data.frame(.basic_stats_round(core, ignore_zeros = ignore_zeros))
  colnames(stats_table) <- if (k == 1) 'Value' else vnames

  combi_plot <- NULL
  if (isTRUE(plot)) {
    panels <- lapply(seq_len(k), function(j) {
      .basic_stats_panel(ts[, j], pstart = pstart, pend = pend, show_label = show_label,
                         label_prefix = label_prefix, show_table = show_table,
                         xpos_label = xpos_label, ypos_label = ypos_label,
                         xpos_table = xpos_table, ypos_table = ypos_table, nbins = nbins,
                         ignore_zeros = ignore_zeros, zero_threshold = zero_threshold)
    })
    if (k == 1) {
      combi_plot <- panels[[1]]
    } else {
      names(panels) <- vnames
      combi_plot <- panels
    }
  }

  list_out <- list(plot = combi_plot, stats_table = stats_table)

  return(list_out)

}


# Column-wise statistics engine. Returns an UNROUNDED matrix with one row per
# statistic (rownames are the stat labels) and one column per input series.
# When ignore_zeros = TRUE the distributional statistics are computed on the
# non-zero values (sub-threshold values masked to NA) and the transition
# statistics are omitted, mirroring the original single-series behaviour.
.basic_stats_core <- function(m, zero_threshold = 0.01, ignore_zeros = FALSE) {
  storage.mode(m) <- 'double'
  n <- nrow(m); k <- ncol(m)

  NumofData         <- rep(n, k)
  NumofMisData      <- matrixStats::colSums2(is.na(m))
  PercOfMissingData <- NumofMisData / NumofData * 100
  Pdr               <- matrixStats::colMeans2(m <= zero_threshold, na.rm = TRUE)

  ms <- m
  if (isTRUE(ignore_zeros)) ms[ms <= zero_threshold] <- NA

  nobs  <- matrixStats::colSums2(!is.na(ms))
  Mean  <- matrixStats::colMeans2(ms, na.rm = TRUE)
  Var   <- matrixStats::colVars(ms,  na.rm = TRUE)
  StDev <- matrixStats::colSds(ms,   na.rm = TRUE)
  Variation <- StDev / Mean
  Min   <- suppressWarnings(matrixStats::colMins(ms, na.rm = TRUE))
  Max   <- suppressWarnings(matrixStats::colMaxs(ms, na.rm = TRUE))

  # Central moments derived from the column-centred matrix (no direct reducer).
  xc <- sweep(ms, 2, Mean, FUN = '-')
  m2 <- matrixStats::colMeans2(xc^2, na.rm = TRUE)
  m3 <- matrixStats::colMeans2(xc^3, na.rm = TRUE)
  m4 <- matrixStats::colMeans2(xc^4, na.rm = TRUE)
  Skewness <- m3 / m2^1.5                       # == moments::skewness (population)
  Kurtosis <- m4 / m2^2                         # == moments::kurtosis (non-excess)
  Mom3     <- matrixStats::colSums2(xc^3, na.rm = TRUE) / (nobs - 1)

  qs <- suppressWarnings(matrixStats::colQuantiles(
    ms, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE, drop = FALSE))
  Q5 <- qs[, 1]; Q25 <- qs[, 2]; Q50 <- qs[, 3]; Q75 <- qs[, 4]; Q95 <- qs[, 5]
  IQR <- abs(Q75 - Q25)

  # L-moments have no matrixStats equivalent -> per-column samlmu (see lmom_stats).
  Lmat <- vapply(seq_len(k), function(j) .basic_stats_lmom(ms[, j]), numeric(7))

  rows <- rbind(
    NumofData = NumofData, NumofMisData = NumofMisData,
    PercOfMissingData = PercOfMissingData,
    Min = Min, Max = Max, Mean = Mean, Var = Var, StDev = StDev,
    Variation = Variation, Mom3 = Mom3, Skewness = Skewness, Kurtosis = Kurtosis,
    Lmean = Lmat[1, ], LScale = Lmat[2, ], L3 = Lmat[3, ], L4 = Lmat[4, ],
    LVariation = Lmat[5, ], LSkewness = Lmat[6, ], LKurtosis = Lmat[7, ],
    Pdr = Pdr, Q5 = Q5, Q25 = Q25, Q50 = Q50, Q75 = Q75, Q95 = Q95, IQR = IQR)

  if (!isTRUE(ignore_zeros)) {
    ct <- .basic_stats_transition(m, zero_threshold)
    rows <- rbind(rows,
      MeanDAfterZero = ct$MeanDAfterZero, VarDAfterZero = ct$VarDAfterZero,
      MeanDBeforeZero = ct$MeanDBeforeZero, VarDBeforeZero = ct$VarDBeforeZero,
      MeanDAfterD = ct$MeanDAfterD, VarDAfterD = ct$VarDAfterD,
      ProbDD = ct$ProbDD, ProbNDND = ct$ProbNDND)
  }

  colnames(rows) <- colnames(m)
  rows
}

# Per-column L-moments: LMean, LScale, L3, L4, L-CV, L-Skewness, L-Kurtosis.
.basic_stats_lmom <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 1) return(rep(NA_real_, 7))
  tryCatch({
    l  <- suppressWarnings(lmom::samlmu(x, nmom = 4, ratios = FALSE, trim = 0))
    lr <- suppressWarnings(lmom::samlmu(x, nmom = 4, ratios = TRUE,  trim = 0))
    c(l[1], l[2], l[3], l[4], l[2] / l[1], lr[3], lr[4])
  }, error = function(e) rep(NA_real_, 7))
}

# Vectorised wet/dry transition statistics across all columns. The lag masks and
# the -99 sentinel reproduce the original single-series logic exactly, including
# NA propagation and the 0/0 = NaN behaviour of the transition probabilities.
.basic_stats_transition <- function(m, thr) {
  n <- nrow(m); k <- ncol(m)
  na_k <- rep(NA_real_, k)
  if (n < 2) {
    return(list(MeanDAfterZero = na_k, VarDAfterZero = na_k,
                MeanDBeforeZero = na_k, VarDBeforeZero = na_k,
                MeanDAfterD = na_k, VarDAfterD = na_k,
                ProbDD = na_k, ProbNDND = na_k))
  }
  Ypos <- m[2:n, , drop = FALSE]
  Ylag <- m[1:(n - 1), , drop = FALSE]

  # Keep the value where the transition fires, NA otherwise. Starting from the
  # numeric value matrix (rather than ifelse, whose result is logical when no
  # element ever fires) keeps the result numeric. `keep %in% TRUE` is TRUE only
  # for genuine TRUEs, so FALSE and NA both map to NA - matching ifelse(., ., NA).
  keep_value <- function(value, keep) {
    value[!(keep %in% TRUE)] <- NA
    value
  }
  zAZ <- keep_value(Ypos, Ypos > thr & Ylag <= thr)  # non-zero after a zero
  zBZ <- keep_value(Ylag, Ypos <= thr & Ylag > thr)  # value preceding a zero
  zDD <- keep_value(Ypos, Ypos > thr & Ylag > thr)   # non-zero after a non-zero

  condDD <- (Ypos > thr) & (Ylag > thr)
  zDDp   <- ifelse(condDD, Ypos, -99)
  condND <- (Ypos == thr) & (Ylag == thr)
  zNDp   <- ifelse(condND, Ypos, -99)

  list(
    MeanDAfterZero  = matrixStats::colMeans2(zAZ, na.rm = TRUE),
    VarDAfterZero   = matrixStats::colVars(zAZ,  na.rm = TRUE),
    MeanDBeforeZero = matrixStats::colMeans2(zBZ, na.rm = TRUE),
    VarDBeforeZero  = matrixStats::colVars(zBZ,  na.rm = TRUE),
    MeanDAfterD     = matrixStats::colMeans2(zDD, na.rm = TRUE),
    VarDAfterD      = matrixStats::colVars(zDD,  na.rm = TRUE),
    ProbDD   = matrixStats::colSums2(zDDp != -99 & !is.na(zDDp)) /
               matrixStats::colSums2(!is.na(zDDp)),
    ProbNDND = matrixStats::colSums2(zNDp != -99 & !is.na(zNDp)) /
               matrixStats::colSums2(!is.na(zNDp)))
}

# Apply the original per-statistic rounding to the unrounded core matrix.
# PercOfMissingData is left unrounded (as in the original).
.basic_stats_round <- function(core, ignore_zeros = FALSE) {
  labs   <- rownames(core)
  digits <- stats::setNames(rep(2L, length(labs)), labs)
  digits[c('NumofData', 'NumofMisData', 'Kurtosis')] <- 0L
  if (isTRUE(ignore_zeros)) {
    digits[c('Q5', 'Q25', 'Q50', 'Q75', 'Q95', 'IQR')] <- 5L
  } else {
    digits[c('MeanDAfterZero', 'VarDAfterZero', 'MeanDBeforeZero', 'VarDBeforeZero',
             'MeanDAfterD', 'VarDAfterD', 'ProbDD', 'ProbNDND')] <- 5L
  }
  out <- core
  for (i in seq_along(labs)) {
    if (labs[i] == 'PercOfMissingData') next
    out[i, ] <- round(core[i, ], digits[[labs[i]]])
  }
  out
}

# Build the timeseries / PDF / ECDF / ACF panel for a single-column series.
.basic_stats_panel <- function(ts, pstart, pend, show_label, label_prefix, show_table,
                               xpos_label, ypos_label, xpos_table, ypos_table, nbins,
                               ignore_zeros, zero_threshold) {
  if (ignore_zeros == TRUE){
    temp <- ts[ts>zero_threshold,]
  }else{
    temp <- ts
  }
  Mean <- round(mean(temp,na.rm=TRUE),2)
  StDev <- round(stats::sd(temp,na.rm=TRUE),2)
  Variation <- round(StDev/Mean,2)
  Mom3 <- round(sum((temp-Mean)^3,na.rm=T)*(1/(sum(!is.na(temp))-1)),6)
  Skewness<-round(Mom3/StDev^3,2)
  Pdr<-round(mean(ts<=zero_threshold,na.rm=T),2)

  stats_basic <- cbind(Metric=c('Mean','Variation','Skewness','Prob. Dry'), Value=c(Mean,Variation,Skewness,Pdr))

  if (!is.na(pstart)){start_date <- pstart}else{start_date <- index(ts)[1]}
  if (!is.na(pend)){end_date <- pend}else{end_date <- index(ts)[nrow(ts)]}
  plot_period <- paste(start_date,end_date,sep = '/')
  ts_plot <- ts[plot_period]

  data_min <- matrixStats::colMins(ts_plot, na.rm = TRUE)
  data_max <- matrixStats::colMaxs(ts_plot, na.rm = TRUE)

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
    geom_histogram(aes(y=after_stat(density)), bins = nbins,     # Histogram with density instead of count on y-axis
                   colour='black', fill='white') +
    geom_density(alpha=.2, fill='#FF6666') + # Overlay with transparent density plot
    labs(x = colnames(ts), y = 'Density') + ggtitle('PDF')


  plot_ecdf <- ggplot(ts_df, aes(X)) + stat_ecdf(geom = 'step') + scale_y_continuous(trans = 'log10', breaks = trans_breaks('log10', function(x) 10^x),
                                                                                     labels = trans_format('log10', math_format(10^.x))) +
    labs(x = colnames(ts), y = 'P(X<x)') + ggtitle('ECDF')

  ts_clean <- ts[!is.na(ts),]
  if (ignore_zeros == TRUE){ts_clean <- ts_clean[which(ts_clean > zero_threshold),]}
  acf_dirty <- stats::acf(ts_clean, lag.max = 10, plot = FALSE)
  acf2df <- data.frame(ACF=acf_dirty$acf, Lag = seq(from = 0, length.out = nrow(acf_dirty$acf)))

  plot_acf <- ggplot(acf2df, aes(x = Lag, y = ACF)) + geom_point() + geom_line() + ggtitle('ACF') #+ ylab(expression(rho_{t,t-1}))

  combi_plot <- (plot_rawts + plot_hist +  plot_layout(widths = c(1.5, 1)))/(plot_ecdf + plot_acf)

  return(combi_plot)
}
