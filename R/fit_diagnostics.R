#' @title fit_diagnostics
#'
#' @description Performs a goodness-of-fit test to diagnose fitted distributions.
#'
#' @param ts A xts object containing the time series data.
#' @param distr The distribution function to be tested.
#' @param params A list object containing the parameters of the fitted distribution.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @return A list containing a Q-Q plot and a P-P plot and a goodness-of-fit table.
#' The GoF metric calculated are the Kramer von-Mises and Kolmogorov-Smirnov. NEEDS EXPANSION
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#'fit <- fitlm_multi(data[,4],candidates = 'gamma3', ignore_zeros = TRUE)
#'
#'fcheck <- fit_diagnostics(data[,4], distr = 'gamma3', params = fit$parameter_list[[1]], ignore_zeros = TRUE)
#'
#' @export
#'
#'



#Function that plots distribution fitting diagnostics, i.e. probability plot, qq and pp plots.
#It requires the names of the quantile and cdf functions and a list of parameters as an argument

fit_diagnostics <- function(ts,dist = 'norm', params,ignore_zeros = FALSE,
                            zero_threshold = 0.01){

  qfunction<-paste0('q',dist)
  pfunction<-paste0('p',dist)
  x <- na.omit(coredata(ts))
  if (ignore_zeros == TRUE){
    if (length(which(x>zero_threshold))>0){
      I <- which(x>zero_threshold)
      pzero <- length(I)/length(x)
      x <- x[I]
    }
  }

  breaks <- c(0,0.25,0.5,0.75,0.9,0.99,0.999,0.9999)
  transform_fun <- function(x){
    param_list_temp1 <- params
    param_list_temp1$p <- x
    do.call(eval(parse(text = qfunction)),param_list_temp1)
  }
  transform_inv <- function(x){
    param_list_temp2 <- params
    x <- ifelse(x<0,0.0001,x)
    #x <- ifelse(is.na(x),0.0001,x)
    param_list_temp2$q <- x
    do.call(eval(parse(text = pfunction)),param_list_temp2)
  }

  fit_trans <- trans_new('trans_fit',transform =  transform_fun,
                         inverse = transform_inv, breaks = breaks, format = format_format(), domain = c(0, 1))

  u_emp <- ppoints(sort(x))
  q_emp <- sort(x)
  emp_data <- data.frame(u_emp = u_emp, q_emp = q_emp)
  u <- seq(0.001,ifelse(max(u_emp)>0.9999,0.9999,max(u_emp)),by = 0.0001)
  param_list_temp3 <- params
  param_list_temp3$p <- u
  q_fitted <- do.call(eval(parse(text = qfunction)),param_list_temp3)
  fit_data <- data.frame(u = u , q_fitted = q_fitted)

  # Probplot <- ggplot() + geom_point(data = emp_data, aes(x=q_emp , y=u_emp, colour = 'Empirical'), size = 2) +
  #   geom_line(data = fit_data, aes(x=q_fitted, y=u, colour = 'Fitted'), size = 1) +
  #   theme(legend.position='bottom') +
  #   labs(colour = 'Legend', x = 'data', y = 'P(x<X)') + scale_colour_manual(values = c('red','black')) +
  #   ggtitle('Probability Plot') +  scale_y_continuous(breaks = breaks, labels = breaks, trans = fit_trans)
  Probplot <- ggplot() + geom_point(data = emp_data, aes(x=q_emp , y=u_emp), shape = 1, size = 2, stroke = 2) +
    geom_line(data = fit_data, aes(x=q_fitted, y=u), size = 1) +
    labs(colour = 'Legend', x = 'data', y = 'P(x<X)') +
    ggtitle('Probability Plot') +  scale_y_continuous(breaks = breaks, labels = breaks, trans = fit_trans)
  #Probplot <- ggplot()
  param_list_temp4 <- params
  param_list_temp4$p <- u_emp
  qq_fitted <- do.call(eval(parse(text = qfunction)),param_list_temp4)

  param_list_temp5 <- params
  param_list_temp5$q <- q_emp
  pp_fitted <- do.call(eval(parse(text = pfunction)),param_list_temp5)

  qq_data <- data.frame(emp = q_emp, fit = qq_fitted)
  qq_line <- data.frame(x = seq(min(x),max(x), by = (max(x)-min(x))/1000), y = seq(min(x),max(x), by = (max(x)-min(x))/1000))

  QQplot <- ggplot() + geom_point(data = qq_data, aes(x=fit, y=emp), shape = 1, size = 2, stroke = 2) +
    geom_line(data = qq_line, aes(x=x, y=y), size = 1) + labs(x = 'Theoretical Quantiles', y = 'Empirical Quantiles') +
    ggtitle('Q-Q Plot')

  pp_data <- data.frame(emp = u_emp, fit = pp_fitted)
  pp_line <- data.frame(x = seq(0,1, by = 0.05), y = seq(0,1, by = 0.05))

  PPplot <- ggplot() + geom_point(data = pp_data, aes(x=fit, y=emp), shape = 1, size = 2, stroke = 2) + xlim(0,1) + ylim(0,1) +
    geom_line(data = pp_line, aes(x=x, y=y), size = 1) + labs(x = 'Theoretical Probabilities', y = 'Empirical Probabilities') +
    ggtitle('P-P Plot')

  diagnostics <- ggpubr::ggarrange(Probplot, QQplot, PPplot, nrow = 1, ncol = 3)

  #MLE <- sum(log(doeval(parse(text = dfunction))(x,par[1],par[2])))
  #KS <- ks.test(x = q_emp, y = qq_fitted)
  CM <- CDFt::CramerVonMisesTwoSamples(q_emp, qq_fitted)
  KS <- CDFt::KolmogorovSmirnov(q_emp, qq_fitted)


  GoF <- list(CramerVonMises = CM,KolmogorovSmirnov = KS)

  out_list <- list(Diagnostic_Plots = diagnostics, GoF = GoF, Probplot = Probplot , QQplot = QQplot, PPplot = PPplot)
  return(out_list)
}
