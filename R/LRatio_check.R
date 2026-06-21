#' Check Sample L-Ratios Against Distribution Support
#'
#' @description Tests whether the sample L-CV and L-skewness of one or more time
#' series lie within the admissible L-ratio space of four 3-parameter
#' distributions: Dagum, generalised gamma, exponentiated Weibull, and Burr type
#' XII. Each distribution's theoretical L-ratio support is represented by a
#' pre-computed polygon. The function identifies the sample points that fall within
#' the theoretical space of each distribution and reports the percentage falling inside each
#' distribution's support. A diagnostic plot overlays the sample L-ratios (red
#' points) on the theoretical L-ratio boundaries, with separate panels for the
#' four 3-parameter distributions and a fifth panel showing the 2-parameter
#' theoretical L-skewness-vs-L-CV curves as a reference. This tool helps
#' pre-screen which distributions are geometrically capable of representing a
#' given dataset before fitting.
#'
#' @param ts_lmoms A data frame or matrix of L-moment statistics as returned by
#'   \code{\link{lmom_stats}}, with rows including L-Skew (row 3) and L-CV
#'   (row 5). Multiple columns are treated as separate series.
#'
#' @return A list with elements \code{distributions} (per-distribution check
#'   results and individual plots) and \code{multi_plots} (a combined
#'   \pkg{patchwork} plot of all panels).
#'
#' @examples
#' # Synthetic daily data
#' set.seed(42)
#' ts <- xts::xts(cbind(series1 = rgamma(3650, shape = 2, scale = 5),
#'                 series2 = rgamma(3650, shape = 3, scale = 7)),
#'           order.by = seq.Date(as.Date("2000-01-01"), by = "day", length.out = 3650))
#' ts_lmoms <- lmom_stats(ts, ignore_zeros = TRUE)
#' lcheck <- LRatio_check(ts_lmoms)
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

    poly_df <- data.frame(x = dist@polygons[[1]]@Polygons[[1]]@coords[,1],
                          y = dist@polygons[[1]]@Polygons[[1]]@coords[,2])

    plot <- ggplot(poly_df, aes(x = x, y = y))+
      ggplot2::geom_polygon(colour='black', fill=NA)+xlab('L-CV') + ggplot2::ylab('L-Skewness')+
      geom_point(data=ltest, aes(CV,skew),color='red') + ylim(-0.5,1) +
      annotate('label', y = -0.1, x = 0.85,label = paste0(sprintf('In = %2.0f',perc_in),'%'))

    list(check=check,plot=plot)

  }


  p_simple <- ggplot(data = LMom_LSpace, aes(x = `L-CV`, y = `L-Skewness`, color = Dist)) + geom_line(linewidth = 1.0) +
    geom_point(data=ltest, aes(CV,skew),color='red') + ggtitle("2-parametric distributions")

  check_space <-  lapply(r_space,FUN = lratio_check, ltest=ltest)
  check_space$Dagum$plot <- check_space$Dagum$plot + ggtitle(names(r_space)[1])
  check_space$GGamma$plot <- check_space$GGamma$plot + ggtitle(names(r_space)[2])
  check_space$ExpWeibull$plot <- check_space$ExpWeibull$plot + ggtitle(names(r_space)[3])
  check_space$BurrXII$plot <- check_space$BurrXII$plot + ggtitle(names(r_space)[4])
  check_space$other$plot <- p_simple

  vis_space_all <-  patchwork::wrap_plots(plotlist  = lapply(check_space,function(x) x$plot))

  list(distributions=check_space, multi_plots = vis_space_all)
}

