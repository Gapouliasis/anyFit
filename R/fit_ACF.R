#' @title fit_ACF
#'
#' @description Function to fit theoretical Autocorrelation Functions (ACFs) to timeseries data.
#' Three functions are available, the Cauchy-type autocorrelation structure,
#' the Hurst - Kolmogorov and the Short Range dependence structure.
#'
#' @param ts A xts object containing the time series data.
#' @param lag_max Maximum lag to use in fitting.
#' @param type A list of character strings, containing the type of ACF to fit.
#'   Options include: \code{'CAS'}, \code{'HK'}, and \code{'SRD'}.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @return A list containing a vector of fitted parameters, a data frame of
#'   fitted ACF values, and a plot of the fitted ACFs.
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#' acfs <- fit_ACF(data[,4], lag = 10, ignore_zeros = TRUE )
#' acfs$ACF_plot
#'
#' @export
#'

fit_ACF <- function(ts,lag_max,type = list('CAS','HK','SRD'),ignore_zeros = FALSE,
                    zero_threshold = 0.01){
  #Omits NANs
  ts <- na.omit(coredata(ts))
  if (ignore_zeros == TRUE){
    ts <- ts[ts > zero_threshold]
  }
  a = stats::acf(ts,lag.max = lag_max, plot = FALSE)
  auto = data.frame(auto=as.numeric(a$acf),
                    lag = a$lag)

  rmse_acf = function(params,auto,lag_max){
    ACF_FUN = paste0(tt,'_ACF')
    param_list = as.list(params)
    param_list$lag_max = lag_max
    names(param_list) = ACF_args
    p_acf <- do.call(eval(parse(text = ACF_FUN)), param_list)
    rmse <- sqrt((sum(p_acf$ACF - auto$auto)^2)/(lag_max + 1))
    return(rmse)
  }

  for (tt in type){
    ACF_FUN = paste0(tt,'_ACF')
    ACF_args <- formalArgs(eval(parse(text = ACF_FUN)))
    nargs = length(ACF_args)
    start = rep(0.5,(nargs-1))
    lbs = rep(0.0001,(nargs-1))
    ubs = rep(Inf,(nargs-1))

    acf_params <- optim(par=start, fn=rmse_acf, auto=auto, lag_max=lag_max,
                        lower=lbs, upper=ubs, method = 'L-BFGS-B')

    param_list = as.list(acf_params$par)
    param_list$lag_max = lag_max
    names(param_list) = ACF_args
    temp <- do.call(eval(parse(text = ACF_FUN)), param_list)
    if (tt == type[1]){
      Iparam = names(param_list) != 'lag_max'
      ACF_params <- list(param_list[Iparam])
      ACF_fitted <- temp
      colnames(ACF_fitted)[2] <- tt
    }else{
      Iparam = names(param_list) != 'lag_max'
      ACF_params <- c(ACF_params, list(param_list[Iparam]))
      temp = as.data.frame(temp[,2])
      names(temp) = tt
      ACF_fitted <- cbind(ACF_fitted, temp)
    }
  }

  names(ACF_params) = type

  long_ACF <- reshape2::melt(ACF_fitted, id = c('lag'))

  ACF_plot <- ggplot() + geom_point(data = auto, aes(x=lag,y=auto),shape = 1, size = 2, stroke = 2) +
    geom_line(data = long_ACF, aes(x=lag, y=value, color = variable),size = 1) +
    ylab('ACF') + labs(color = 'Legend') + scale_color_brewer(palette='Set1')

  return(list('ACF_params' = ACF_params, 'ACF_fitted' = ACF_fitted
              , 'ACF_plot' = ACF_plot))
}


