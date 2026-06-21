#' @title Multi-Distribution Fitting
#'
#' @description Fits a list of candidate distributions to a single time series
#'   via the L-moments method. Each candidate is dispatched to its corresponding
#'   \code{fitlm_*} function, which returns the fitted parameters, theoretical
#'   L-moments, sample L-moments, and six goodness-of-fit metrics (MLE, CM, KS,
#'   MSEquant, DiffOfMax, MeanDiffOf10Max). When \code{diagnostic_plots = TRUE}, Q-Q and P-P
#'   comparison plots are produced with all candidates overlaid on a single panel,
#'   using colour-coded points and a 1:1 reference line. An empty series or one with all values
#'   below the zero threshold after filtering returns an \code{NA} placeholder
#'   for every candidate.
#'
#' @param ts An xts object containing the time series data.
#' @param candidates A character vector of distribution names to fit.
#' @param ignore_zeros A logical value, if \code{TRUE} zeros will be ignored.
#'   Default is \code{FALSE}.
#' @param zero_threshold The threshold below which values are considered zero.
#'   Default is 0.01.
#' @param diagnostic_plots A logical value controlling whether Q-Q and P-P
#'   diagnostic plots are produced. Default is \code{TRUE}.
#' @param order Optional named list mapping a candidate name to the vector of
#'   L-moment orders matched by its optimiser, e.g.
#'   \code{list(gengamma = 1:5, expweibull = 1:3)}. Only the numerically-fitted
#'   distributions (\code{gengamma}, \code{gengamma_loc}, \code{burr},
#'   \code{dagum}, \code{expweibull}) accept it; entries for other candidates are
#'   ignored. Candidates not named keep their own default. Default \code{NULL}
#'   (every distribution uses its default).
#'
#' @return A list with components \code{parameter_list} (a named list of
#'   per-candidate fit results, each containing \code{Distribution},
#'   \code{Param}, \code{TheorLMom}, \code{DataLMom}, and \code{GoF}),
#'   \code{GoF_summary} (a data frame of GoF metrics with candidates as columns),
#'   and, when \code{diagnostic_plots = TRUE}, \code{diagnostics} (a combined
#'   Q-Q and P-P panel), \code{QQplot}, and \code{PPplot}.
#'
#' @examples
#' # Daily precipitation-like data: gamma-distributed with zeros
#' x <- xts::xts(rgamma(365, shape = 0.8, scale = 3),
#'          order.by = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 365))
#' x[sample(1:365, 100)] <- 0
#'
#' candidates <- c('exp', 'gamma3', 'weibull')
#' fits <- fitlm_multi(x, candidates = candidates, ignore_zeros = TRUE)
#' fits$diagnostics
#' fits$GoF_summary
#'
#' @export
#'

fitlm_multi <- function(ts,candidates,ignore_zeros = FALSE, zero_threshold = 0.01, diagnostic_plots = TRUE, order = NULL){
  x <- stats::na.omit(coredata(ts))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }
  if (length(x) == 0) {
    empty_fit <- list(
      Distribution = list(FXs = NA_character_),
      Param = NA_real_,
      TheorLMom = rep(NA_real_, 4),
      DataLMom = rep(NA_real_, 4),
      GoF = list(MLE = NA_real_, CM = NA_real_, KS = NA_real_,
                 MSEquant = NA_real_, DiffOfMax = NA_real_, MeanDiffOf10Max = NA_real_)
    )
    params_list <- stats::setNames(
      lapply(candidates, function(c) empty_fit),
      as.character(candidates)
    )
    GoF_df <- as.data.frame(matrix(NA_real_, nrow = 6, ncol = length(candidates)))
    rownames(GoF_df) <- c("MLE", "CM", "KS", "MSEquant", "DiffOfMax", "MeanDiffOf10Max")
    colnames(GoF_df) <- as.character(candidates)
    if (diagnostic_plots) {
      return(list(parameter_list = params_list, GoF_summary = GoF_df,
                  diagnostics = NULL, QQplot = NULL, PPplot = NULL))
    } else {
      return(list(parameter_list = params_list, GoF_summary = GoF_df))
    }
  }
  if (diagnostic_plots){
    u_emp <- ppoints(sort(x))
    q_emp <- sort(x)
    emp_data <- data.frame(u_emp = u_emp, q_emp = q_emp)
  }

  for (candidate in candidates){
    fit_function <- match.fun(paste0('fitlm_',candidate))
    fit_args <- list(x = x, ignore_zeros = FALSE, zero_threshold = 0.01)
    if (!is.null(order[[candidate]]) && "order" %in% names(formals(fit_function))) {
      fit_args$order <- order[[candidate]]
    }
    params <- do.call(fit_function, fit_args)

    if (candidate == candidates[[1]]){
      params_list <- list(params)
      GoF <- as.data.frame(unlist(params$GoF))
      colnames(GoF) <- candidate
    }else{
      params_list <- c(params_list,list(params))
      temp <- as.data.frame(unlist(params$GoF))
      colnames(temp) <- candidate
      GoF <- cbind(GoF, temp)
    }

    if (diagnostic_plots){
      param_list_temp1 <- params$Param
      param_list_temp1$p <- u_emp
      qq_fitted <- do.call(match.fun(paste0('q',candidate)),param_list_temp1)

      param_list_temp2 <- params$Param
      param_list_temp2$q <- q_emp
      pp_fitted <- do.call(match.fun(paste0('p',candidate)),param_list_temp2)

      if (candidate == candidates[[1]]){
        qq_data <- data.frame(emp = q_emp, fit = qq_fitted, FX = candidate)
        pp_data <- data.frame(emp = u_emp, fit = pp_fitted, FX = candidate)
      }else{
        temp1 <- data.frame(emp = q_emp, fit = qq_fitted, FX = candidate)
        qq_data <- rbind(qq_data, temp1)
        temp2 <- data.frame(emp = u_emp, fit = pp_fitted, FX = candidate)
        pp_data <- rbind(pp_data, temp2)
      }
    }

  }

  names(params_list) <- unlist(candidates)

  if (diagnostic_plots){
    qq_line <- data.frame(x = seq(min(x),max(x), by = (max(x)-min(x))/1000), y = seq(min(x),max(x), by = (max(x)-min(x))/1000))
    pp_line <- data.frame(x = seq(0,1, by = 0.05), y = seq(0,1, by = 0.05))

    QQplot <- ggplot() + geom_point(data = qq_data, aes(x=fit, y=emp, color = FX), shape = 1, size = 1.5, stroke = 1.5) +
      geom_line(data = qq_line, aes(x=x, y=y), linewidth = 1) + labs(x = 'Theoretical Quantiles', y = 'Empirical Quantiles') +
      ggtitle('Q-Q Plot') + scale_color_brewer(palette='Set1')

    PPplot <- ggplot() + geom_point(data = pp_data, aes(x=fit, y=emp, color = FX), shape = 1, size = 1.5, stroke = 1.5) + xlim(0,1) + ylim(0,1) +
      geom_line(data = pp_line, aes(x=x, y=y), linewidth = 1) + labs(x = 'Theoretical Probabilities', y = 'Empirical Probabilities') +
      ggtitle('P-P Plot') + scale_color_brewer(palette='Set1')

    combined <-  patchwork::wrap_plots(QQplot,PPplot, nrow = 1, ncol = 2)

    list_out <- list('parameter_list' = params_list, 'GoF_summary' = GoF,
                     'diagnostics' = combined, 'QQplot' = QQplot, 'PPplot' = PPplot)
  }else{
    list_out <- list('parameter_list' = params_list, 'GoF_summary' = GoF)
  }

  return(list_out)
}
