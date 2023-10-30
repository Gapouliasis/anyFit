library(anyFit)
filename = "C:/Users/Admin/Downloads/rr_ens_mean_0.25deg_reg_2011-2022_v27.0e.nc"

data = nc2xts(filename = filename, varname = "rr", country = "Belgium")
data_raster = data$raster

annual_data = period_apply_nc(data_raster, period = "years", FUN = "sum")

annual_plot = nc_ggplot(annual_data,title = TRUE,viridis.option = "turbo",
                        legend.title = "Annual Precipitation (mm)", common.legend = TRUE, legend = "bottom")

stats = basic_stats_nc(annual_data, ignore_zeros = TRUE)

stats_plot = nc_ggplot(stats[[4:15]])

fits_all = fitlm_nc(data_raster, ignore_zeros = TRUE, candidates = "gamma")

params_plot = nc_ggplot(fits_all$raster_params, viridis.option = "turbo")

Lmom_plot = nc_ggplot(fits_all$raster_TheorLMom, viridis.option = "inferno")

nc_data = ncdf4::nc_open(filename = filename)

y_coords = ncdf4::ncvar_get(nc_data, "latitude")
x_coords = ncdf4::ncvar_get(nc_data, "longitude")

coords = data.frame(x = x_coords[c(220)],y = y_coords[c(150)])

rnd_xts = nc2xts_nn(filename = filename, varname = "rr", coords = coords)

params = fitlm_gengamma(data2[,4], ignore_zeros = TRUE)
params = fitlm_burr(rnd_xts, ignore_zeros = TRUE)
fit_check <- fit_diagnostics(rnd_xts, dist = 'burr',
                             params = params$Param$para, ignore_zeros = TRUE)

params = fitlm_dagum(rnd_xts, ignore_zeros = TRUE)
fit_check <- fit_diagnostics(rnd_xts, dist = 'dagum',
                             params = params$Param, ignore_zeros = TRUE)

params = fitlm_gengamma(rnd_xts, ignore_zeros = TRUE)
fit_check <- fit_diagnostics(rnd_xts, dist = 'gengamma',
                             params = params$Param, ignore_zeros = TRUE)

params = fitlm_gengamma_loc(rnd_xts, location = 0, ignore_zeros = TRUE)
