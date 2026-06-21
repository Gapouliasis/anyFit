#' @title fit_ACF
#'
#' @description Fits one or more theoretical autocorrelation functions to
#'   the sample ACF of an xts series. Three parametric models are
#'   available: the Cauchy-type autocorrelation structure (CAS, a
#'   two-parameter model with algebraic decay), the Hurst-Kolmogorov
#'   process (HK, a one-parameter fractional Gaussian noise model), and
#'   the short-range dependence model (SRD, a one-parameter exponential
#'   decay). Parameters are estimated by minimising the root-mean-square
#'   error between the theoretical and sample ACFs up to
#'   \code{lag_max}, using the L-BFGS-B optimiser with positivity
#'   constraints on all parameters. The function returns fitted parameter
#'   vectors for each model type, a data frame of theoretical ACF values at
#'   each lag, and a \code{ggplot} overlay of the fitted curves on the
#'   empirical autocorrelation.
#'
#' @param ts An xts object containing the time series data.
#' @param lag_max Integer. Maximum lag to use in fitting.
#' @param type Character vector of ACF model types. Options are
#'   \code{"CAS"}, \code{"HK"}, and \code{"SRD"}. Default all three.
#' @param ignore_zeros Logical. If \code{TRUE}, values below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric. Threshold below which values are treated
#'   as zero. Default \code{0.01}.
#'
#' @return A list with elements \code{ACF_params} (named list of fitted
#'   parameter vectors per model), \code{ACF_fitted} (data frame of
#'   theoretical ACF values with a \code{lag} column), and
#'   \code{ACF_plot} (\code{ggplot} of sample ACF with fitted curves
#'   overlaid).
#'
#' @examples
#' ts <- xts::xts(rnorm(365), order.by = as.Date("2020-01-01") + 0:364)
#' acfs <- fit_ACF(ts, lag_max = 10, type = c("CAS", "HK", "SRD"))
#' acfs$ACF_plot
#'
#' @importFrom methods formalArgs
#' @export
#'

fit_ACF <- function(ts,lag_max,type = list('CAS','HK','SRD'),ignore_zeros = FALSE,
                    zero_threshold = 0.01){
  #Omits NANs
  ts <- stats::na.omit(coredata(ts))
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

    acf_params <- stats::optim(par=start, fn=rmse_acf, auto=auto, lag_max=lag_max,
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
    geom_line(data = long_ACF, aes(x=lag, y=value, color = variable),linewidth = 1) +
    ggplot2::ylab('ACF') + labs(color = 'Legend') + scale_color_brewer(palette='Set1')

  return(list('ACF_params' = ACF_params, 'ACF_fitted' = ACF_fitted
              , 'ACF_plot' = ACF_plot))
}

