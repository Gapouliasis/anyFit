i_months <- unique(month(ts))
month_name <- rep(0,length(i_months))
monthly_plots = list()
for (i in i_months){
  month_name[i] <- month.name[i]
  I <- which(month(ts) == i)
  monthly_ts <- ts[I]
  monthly_ts_long = reshape2::melt(coredata(monthly_ts))
  
  monthly_ecdf = apply(coredata(monthly_ts), function(x){rank(coredata(x))}, MARGIN = 2)/nrow(monthly_ts)
  monthly_ecdf_long = reshape2::melt(monthly_ecdf)
  monthly_ecdf_long$q = monthly_ts_long$value
  
  monthly_ecdf_plot = ggplot(monthly_ecdf_long, aes(x = q, y = value, color = Var2)) + 
    geom_point(shape = 1, size = 1.5, stroke = 1.5) + ggtitle(month.name[i]) +
    scale_color_brewer(palette='Set1') + 
    labs(x = 'Empirical Quantile', y = 'Empirical Probability')
    theme(legend.position = 'bottom')
  monthly_plots = c(monthly_plots, list(monthly_ecdf_plot))
}

ggpubr::ggarrange(plotlist = monthly_plots,
                  nrow = 4, ncol = 3,
                  common.legend = TRUE, legend = 'bottom')
