#' Basic statistics and diagnostic panel for time series data
#'
#' @description
#' Computes a comprehensive set of summary statistics for one or more time
#' series stored as xts columns and optionally produces a four-panel diagnostic
#' plot (time series, probability density function, empirical cumulative
#' distribution function, and autocorrelation function). The statistics include
#' sample moments (mean, variance, skewness, kurtosis), L-moments, quantiles,
#' probability dry, and wet/dry transition statistics. When \code{ignore_zeros}
#' is \code{TRUE}, distributional statistics are computed only on values
#' exceeding \code{zero_threshold}, and transition statistics are omitted.
#' Multi-column input is supported: the statistics table then contains one
#' column per input series (named after the xts columns), and the plot output
#' is a named list of one panel per series. The internal computation is
#' vectorised across columns via \code{matrixStats} for computational efficiency.
#'
#' @param ts An xts object. If multi-column, each column is treated as an
#'   independent series.
#' @param pstart Plot start date (character or date). If \code{NA}, the first
#'   date of the series is used.
#' @param pend Plot end date (character or date). If \code{NA}, the last date
#'   of the series is used.
#' @param show_label Logical. If \code{TRUE}, a label with the series name is
#'   plotted on the timeseries panel. Default is \code{TRUE}.
#' @param label_prefix Character prefix for the series label, e.g.
#'   \code{"Station"}. Default is \code{"Station"}.
#' @param show_table Logical. If \code{TRUE}, a table of key statistics (mean,
#'   coefficient of variation, skewness, probability dry) is overlaid on the
#'   timeseries panel. Default is \code{TRUE}.
#' @param xpos_label Horizontal position of the series label, as a fraction of
#'   the plot width (0–1). Default is \code{0.1}.
#' @param ypos_label Vertical position of the series label, as a fraction of
#'   the data range (0–1). Default is \code{0.95}.
#' @param xpos_table Horizontal position of the statistics table, as a fraction
#'   of the plot width (0–1). Default is \code{0.1}.
#' @param ypos_table Vertical position of the statistics table, as a fraction
#'   of the data range (0–1). Default is \code{0.15}.
#' @param nbins Number of bins for the PDF histogram. Default is \code{30}.
#' @param ignore_zeros Logical. If \code{TRUE}, values at or below
#'   \code{zero_threshold} are excluded from distributional statistics and
#'   transition statistics are omitted. Default is \code{FALSE}.
#' @param zero_threshold Numeric threshold below which values are treated as
#'   zero. Default is \code{0.01}.
#' @param plot Logical. If \code{TRUE}, the diagnostic panel is built. Default
#'   is \code{FALSE}. For multi-column input one panel is produced per column.
#'
#' @return A list with two elements: \code{stats_table}, a data frame of
#'   statistics (rows are the statistics, columns are the input series — named
#'   after the input columns when more than one is supplied, otherwise
#'   \code{"Value"}); and \code{plot}, which is \code{NULL} when
#'   \code{plot = FALSE}, a single combined ggplot for a single-column input,
#'   or a named list of one combined ggplot per column for a multi-column input.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density labs ggtitle stat_ecdf scale_y_continuous geom_point geom_line annotate annotation_custom xlab theme element_text after_stat
#' @importFrom patchwork plot_layout
#' @importFrom scales trans_breaks trans_format math_format
#' @importFrom zoo autoplot.zoo coredata index
#' @importFrom matrixStats colMins colMaxs colMeans2 colVars colSds colSums2 colQuantiles
#'
#' @examples
#' # Synthetic daily precipitation series (two stations)
#' set.seed(123)
#' n <- 1000
#' dates <- seq(as.Date("2000-01-01"), by = "day", length.out = n)
#' prec <- xts::xts(cbind(
#'   StationA = rgamma(n, shape = 0.6, scale = 5),
#'   StationB = rgamma(n, shape = 0.7, scale = 4)),
#'   order.by = dates)
#'
#' # Single-column: statistics table only
#' bs <- basic_stats(prec[, 1], ignore_zeros = TRUE)
#' bs$stats_table
#'
#' # Single-column: with diagnostic panel
#' bs <- basic_stats(prec[, 1], pstart = "2001", pend = "2002",
#'                    show_table = FALSE, show_label = FALSE,
#'                    ignore_zeros = TRUE, plot = TRUE)
#' bs$plot
#'
#' # Multi-column: one column of statistics per series, named plot list
#' wbs <- basic_stats(prec, plot = TRUE)
#' wbs$stats_table
#' wbs$plot[[1]]
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


#' Column-wise statistics engine
#'
#' @description
#' Computes an unrounded statistics matrix with one row per statistic and one
#' column per input series. When \code{ignore_zeros = TRUE}, distributional
#' statistics are computed on the non-zero values only (sub-threshold values
#' masked to \code{NA}) and transition statistics are omitted.
#'
#' @param m A numeric matrix (columns = series, rows = observations).
#' @param zero_threshold Numeric threshold. Default is \code{0.01}.
#' @param ignore_zeros Logical. Default is \code{FALSE}.
#'
#' @return A numeric matrix with rownames identifying each statistic and
#'   colnames inherited from \code{m}.
#' @keywords internal
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

#' Per-column L-moments
#'
#' @description
#' Computes L-mean, L-scale, L3, L4, L-CV, L-skewness, and L-kurtosis for a
#' single numeric vector via \code{\link[lmom]{samlmu}}.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector of length 7: LMean, LScale, L3, L4, L-CV,
#'   L-Skewness, L-Kurtosis. Returns \code{NA} for all elements on error or
#'   if \code{x} is empty.
#' @keywords internal
.basic_stats_lmom <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 1) return(rep(NA_real_, 7))
  tryCatch({
    l  <- suppressWarnings(lmom::samlmu(x, nmom = 4, ratios = FALSE, trim = 0))
    lr <- suppressWarnings(lmom::samlmu(x, nmom = 4, ratios = TRUE,  trim = 0))
    c(l[1], l[2], l[3], l[4], l[2] / l[1], lr[3], lr[4])
  }, error = function(e) rep(NA_real_, 7))
}

#' Vectorised wet/dry transition statistics
#'
#' @description
#' Computes per-column transition statistics (mean and variance after/before
#' zero, mean and variance after non-zero, and transition probabilities) using
#' lag masks and a sentinel-based counting scheme.
#'
#' @param m A numeric matrix (columns = series, rows = observations).
#' @param thr Numeric threshold defining the zero/non-zero boundary.
#'
#' @return A list with eight named numeric vectors, each of length
#'   \code{ncol(m)}.
#' @keywords internal
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

#' Apply per-statistic rounding rules
#'
#' @description
#' Rounds each row of the core statistics matrix according to a fixed
#' per-statistic precision table. \code{PercOfMissingData} is left unrounded.
#'
#' @param core A numeric matrix with rownames identifying each statistic.
#' @param ignore_zeros Logical. When \code{TRUE}, quantile rows are rounded to
#'   5 decimal places; otherwise transition statistic rows are rounded to 5
#'   decimal places.
#'
#' @return The rounded numeric matrix (same dimensions as \code{core}).
#' @keywords internal
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

#' Build the diagnostic panel for a single series
#'
#' @description
#' Constructs a four-panel \code{patchwork} plot (raw time series, PDF
#' histogram with density overlay, empirical CDF on log-probability scale, and
#' sample ACF) for a single-column xts series.
#'
#' @param ts A single-column xts object.
#' @param pstart Plot start date (\code{NA} for series start).
#' @param pend Plot end date (\code{NA} for series end).
#' @param show_label Logical. Overlay series label.
#' @param label_prefix Character prefix for the label.
#' @param show_table Logical. Overlay statistics table.
#' @param xpos_label Horizontal label position (0–1 fraction).
#' @param ypos_label Vertical label position (0–1 fraction).
#' @param xpos_table Horizontal table position (0–1 fraction).
#' @param ypos_table Vertical table position (0–1 fraction).
#' @param nbins Number of PDF histogram bins.
#' @param ignore_zeros Logical. Exclude sub-threshold values.
#' @param zero_threshold Numeric zero threshold.
#'
#' @return A \code{patchwork} object combining four ggplot panels.
#' @keywords internal
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
