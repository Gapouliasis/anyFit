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

x_coords = ncdf4::ncvar_get(nc_data, "latitude")
pap
#names(nc_data$var)

nc_brick = raster::brick(filename, varname = varname, level = 1)
