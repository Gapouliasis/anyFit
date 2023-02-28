#' @title correl_plots
#'
#' @description Function that calculates correlation plots for two timeseries.
#'Scatter plot with Pearson correlation coefficient
#'Empirical copula plot with Spearman correlation coefficient
#'Standardized or Normal plot with Pearson semi-correlations
#'Handles missing values (NA) with casewise deletion
#'Can check and retain only common dates. Useful for cross-correlation between different stations and variables.
#'
#'
#' @param x,y A xts objects containing the time series data.
#' @param check_common logical value, if TRUE the function will test for common dates and excludes the rest.
#' Default is TRUE.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @return A list with the combined plots and the individual plots for greater flexibility.
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#' correls <- correl_plots(data[,1],data[,4], check_common = TRUE, ignore_zeros = TRUE)
#'
#' correls$combined
#'
#'
#' @export


correl_plots <- function(x,y, check_common = TRUE,ignore_zeros = FALSE,
                         zero_threshold = 0.01){
  if (check_common == FALSE){
    tryCatch({
      if (ignore_zeros == TRUE){
        x <- x[x > zero_threshold,]
        y <- y[y > zero_threshold,]
      }
      data <- data.frame("X" = coredata(x), "Y" = coredata(y))
    },error = function(e){
      message("Maybe the timeseries have different indexes (dates). Try setting check_common = TRUE")
      message("!Original Error!")
      message(e)
    })

  }else{
    if (ignore_zeros == TRUE){
      x <- x[x > zero_threshold,]
      y <- y[y > zero_threshold,]
    }
    cdata <- merge.xts(x,y, join = "inner")
    data <- as.data.frame(coredata(cdata))
    #code that does the same. More generic
    # a <- data[,1]
    # a <- a[!is.na(a),]
    # b <- data[,3]
    # b <- b[!is.na(b),]
    #
    # mylist <- list(a = time(a), b = time(b))
    #
    # common_values = Reduce(intersect, mylist)
    # ID <- lapply(mylist, function(x) which(x %in% common_values))
    #
    # acommon <- a[ID$a,]
    # bcommon <- b[ID$b,]
    # common <- merge.xts(acommon,bcommon)
  }

  colnames(data) <- c("X","Y")
  rho_pearson <- round(stats::cor(data$X, data$Y, method = "pearson", use = "complete.obs"), digits = 2)

  pso_data <- as.data.frame(copula::pobs(data))
  rho <- round(stats::cor(pso_data$X, pso_data$Y, method = "spearman", use = "complete.obs"), digits = 2)

  x_norm <- qnorm(pso_data$X)
  y_norm <- qnorm(pso_data$Y)

  nrho1 <- round(stats::cor(x_norm[x_norm>=0 & y_norm >=0], y_norm[x_norm>=0 & y_norm >=0], method = "pearson", use = "complete.obs"), digits = 2)
  nrho2 <- round(stats::cor(x_norm[x_norm>=0 & y_norm <=0], y_norm[x_norm>=0 & y_norm <=0], method = "pearson", use = "complete.obs"), digits = 2)
  nrho3 <- round(stats::cor(x_norm[x_norm<=0 & y_norm <=0], y_norm[x_norm<=0 & y_norm <=0], method = "pearson", use = "complete.obs"), digits = 2)
  nrho4 <- round(stats::cor(x_norm[x_norm<=0 & y_norm >=0], y_norm[x_norm<=0 & y_norm >=0], method = "pearson", use = "complete.obs"), digits = 2)

  correl_table <- rbind(rho_pearson, rho, nrho1, nrho2, nrho3, nrho4)
  rownames(correl_table) <- c("Pearson", "Spearman", "Semi_1", "Semi_2", "Semi_3", "Semi_4")

  norm_data <- data.frame("X" = x_norm, "Y" = y_norm)

  dots <- data.frame("X" = seq(from=4 , to = -4 , by = -0.1), "Y" = rep(0,length(seq(from=4 , to = -4 , by = -0.1))))
  label_x <- (max(data$X, na.rm = TRUE) -min(data$X, na.rm = TRUE))*0.15 + min(data$X, na.rm = TRUE)
  label_y <- (max(data$Y, na.rm = TRUE) -min(data$Y, na.rm = TRUE))*1.05 + min(data$Y, na.rm = TRUE)

  plot11 <- ggplot(data = data, aes(x=X, y=Y)) + geom_point( size = 2) +
    scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
    annotate("label", x = label_x, y = label_y, label = sprintf("rho==%1.2f",rho_pearson), parse = TRUE) +
    labs(color = 'Case', y = colnames(y), x = colnames(x)) + ggtitle("Scatter Plot")

  plot12 <- ggplot(data = pso_data, aes(x=X, y=Y)) + geom_point( size = 2) +
    scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
    annotate("label", x = 0.15, y = 1, label = sprintf("r==%1.2f",rho), parse = TRUE) +
    labs(color = 'Case', y = colnames(y), x = colnames(x)) + ggtitle("Copula Plot") + xlim(0, 1)

  plot13 <- ggplot(data = norm_data, aes(x=X, y=Y)) + geom_point( size = 2) +
    geom_line(data = dots, aes(x=X,y=Y), linetype = "dashed") +
    geom_line(data = dots, aes(x=Y,y=X), linetype = "dashed") +
    scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
    annotate("label", x = 2.5, y = 3.5, label = sprintf("rho==%1.2f",nrho1), parse = TRUE) +
    annotate("label", x = 2.5, y = -3.5, label = sprintf("rho==%1.2f",nrho2), parse = TRUE) +
    annotate("label", x = -2.5, y = -3.5, label = sprintf("rho==%1.2f",nrho3), parse = TRUE) +
    annotate("label", x = -2.5, y = 3.5, label = sprintf("rho==%1.2f",nrho4), parse = TRUE) +
    labs(color = 'Case', y = colnames(y), x = colnames(x)) + ggtitle("Standard Normal Plot") + xlim(-3.5, 3.5) + ylim(-3.5,3.5)

  f <- ggpubr::ggarrange(plot11, plot12, plot13, ncol = 3 , nrow = 1)
  output_list <- list(combined = f, scatter_plot = plot11, copula_plot = plot12, normal_plot = plot13, correl_table = correl_table)
  return(output_list)
}




