#' @title fitlm_multi
#'
#' @description This function fits a list of candidate distributions using the L-moments method
#' to a timeseries in xts format. Additionally to the list of the fitted parameters,
#'  goodness-of-fit metric, PP and QQ plots.
#'
#' @param ts A xts object containing the time series data.
#' @param candidates A list of distribution to fit.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @return A list containing the fitted parameters, Goodness-of-Fit Summary, and diagnostic plots.
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#'candidates <- list('exp','expweibull', 'gamma3')
#'fits <- fitlm_multi(data['2010',4],candidates = candidates, ignore_zeros = TRUE)
#'
#' fits$diagnostics
#'
#' @export
#'

fitlm_multi <- function(ts,candidates,ignore_zeros = FALSE, zero_threshold = 0.01){
  x <- na.omit(coredata(ts))
  if (ignore_zeros == TRUE){
    x <- x[x > zero_threshold,]
  }
  u_emp <- ppoints(sort(x))
  q_emp <- sort(x)
  emp_data <- data.frame(u_emp = u_emp, q_emp = q_emp)

  for (candidate in candidates){
    fit_function <- match.fun(paste0('fitlm_',candidate))
    fit_args <- list(x = x, ignore_zeros = FALSE, zero_threshold = 0.01)
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

    qfunction <- paste0('q',candidate)
    pfunction <- paste0('p',candidate)

    param_list_temp1 <- params$Param
    param_list_temp1$p <- u_emp
    qq_fitted <- do.call(eval(parse(text = qfunction)),param_list_temp1)

    param_list_temp2 <- params$Param
    param_list_temp2$q <- q_emp
    pp_fitted <- do.call(eval(parse(text = pfunction)),param_list_temp2)

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


  qq_line <- data.frame(x = seq(min(x),max(x), by = (max(x)-min(x))/1000), y = seq(min(x),max(x), by = (max(x)-min(x))/1000))
  pp_line <- data.frame(x = seq(0,1, by = 0.05), y = seq(0,1, by = 0.05))

  QQplot <- ggplot() + geom_point(data = qq_data, aes(x=fit, y=emp, color = FX), shape = 1, size = 1.5, stroke = 1.5) +
    geom_line(data = qq_line, aes(x=x, y=y), size = 1) + labs(x = 'Theoretical Quantiles', y = 'Empirical Quantiles') +
    ggtitle('Q-Q Plot') + scale_color_brewer(palette='Set1')

  PPplot <- ggplot() + geom_point(data = pp_data, aes(x=fit, y=emp, color = FX), shape = 1, size = 1.5, stroke = 1.5) + xlim(0,1) + ylim(0,1) +
    geom_line(data = pp_line, aes(x=x, y=y), size = 1) + labs(x = 'Theoretical Probabilities', y = 'Empirical Probabilities') +
    ggtitle('P-P Plot') + scale_color_brewer(palette='Set1')

  combined <-  patchwork::wrap_plots(QQplot,PPplot, nrow = 1, ncol = 2)

  names(params_list) <- unlist(candidates)

  return(list('parameter_list' = params_list, 'GoF_summary' = GoF,
              'diagnostics' = combined, 'QQplot' = QQplot, 'PPplot' = PPplot))
}
