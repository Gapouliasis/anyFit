#' @title LRatio_check
#'
#' @description This function checks if the L-ratios of a given time series
#' are within the acceptable parameter space for four different distributions.
#' These are the Dagum, GGamma, ExpWeibull, and BurrXII. TO BE EXPANDED
#'
#' @param ts_lmoms A vector or a matrix of L-moments of the time series of interest.
#' Can test multiple timeseries as columns of a matrix.
#'
#' @return A list containing a boolean (TRUE/FALSE) vector with the results of the test
#' and a ggplot object with the visual representation of the test.
#'
#' @examples
#'
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#' ts_lmoms <- lmom_stats(data, ignore_zeros = TRUE)
#'
#' lcheck <- LRatio_check(ts_lmoms)
#'
#' lcheck$multi_plots
#'
#' @export

LRatio_check <- function(ts_lmoms){
  # Load ratio space -------------------------------------------------------
  # r_space <- list(Dagum=readRDS('Lratio_shp_Dagum.rds'),
  #                 GGamma=readRDS('Lratio_shp_GG.rds'),
  #                 ExpWeibull=readRDS('Lratio_shp_ExpWei.rds'),
  #                 BurrXII=readRDS('Lratio_shp_burrXII.rds'))

  r_space <- list(Dagum=Dagum_LSpace,
                  GGamma=GGamma_LSpace,
                  ExpWeibull=ExpWeibull_LSpace,
                  BurrXII=BurrXII_LSpace)

  # Read L-moments ----------------------------------------------------------
  ltest <- data.frame(CV = as.numeric(ts_lmoms[5,]),
                      skew = as.numeric(ts_lmoms[3,]))

  # Check if the rations <<live>> inside the parameter space and Visualize ----------------
  lratio_check <- function(ltest,dist, dist_name){

    check <- sp::point.in.polygon(point.x = ltest$CV,
                                  point.y = ltest$skew,
                                  pol.x = dist@polygons[[1]]@Polygons[[1]]@coords[,1],
                                  pol.y = dist@polygons[[1]]@Polygons[[1]]@coords[,2])

    perc_in <- round(sum(check)/length(check),2)*100

    plot <- ggplot(dist, aes(x = dist@polygons[[1]]@Polygons[[1]]@coords[,1],
                             y = dist@polygons[[1]]@Polygons[[1]]@coords[,2]))+
      geom_polygon(colour='black', fill=NA)+xlab('L-CV') + ylab('L-Skewness')+
      geom_point(data=ltest, aes(CV,skew),color='red') + ylim(-0.5,1) +
      annotate('label', y = -0.1, x = 0.85,label = paste0(sprintf('In = %2.0f',perc_in),'%'))

    list(check=check,plot=plot)

  }

  check_space <-  lapply(r_space,FUN = lratio_check, ltest=ltest)
  check_space$Dagum$plot <- check_space$Dagum$plot + ggtitle(names(r_space)[1])
  check_space$GGamma$plot <- check_space$GGamma$plot + ggtitle(names(r_space)[2])
  check_space$ExpWeibull$plot <- check_space$ExpWeibull$plot + ggtitle(names(r_space)[3])
  check_space$BurrXII$plot <- check_space$BurrXII$plot + ggtitle(names(r_space)[4])

  vis_space_all <-  patchwork::wrap_plots(plotlist  = sapply(check_space,function(x) x[2]))

  list(test=check_space, multi_plots = vis_space_all)


}


