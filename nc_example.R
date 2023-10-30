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

x <- rnd_xts
ignore_zeros = TRUE
zero_threshold = 0.01
x <- na.omit(coredata(x))
if (ignore_zeros == TRUE){
  NZ=x[x>zero_threshold,]
  PW=length(NZ)/length(x)
  PD=1-PW
}else{
  PW <- 1
  PD <- 0
}

lm=lmom::samlmu(NZ)[1:4]
pfunction=pburr
qfunction=qburr
sample_LM = Lmoments::Lcoefs(NZ, rmax = 4, na.rm = FALSE, trim = c(0, 0))
loc = sample_LM[1]
sample_LM[1] = 0
NZ = NZ - loc
#data_u = apply(data[,c(4,5,6,7)], MARGIN = 2, FUN = function(x){x/max(x)})
I = which.min(sqrt((GenGamma_InitValues[,'L1']-loc)^2+(GenGamma_InitValues[,'L2']-sample_LM[,'L2'])^2+
                     (GenGamma_InitValues[,'tau3']-sample_LM[,'tau3'])^2+
                     (GenGamma_InitValues[,'tau4']-sample_LM[,'tau4'])^2))

I = RANN::nn2(Burr_InitValues[,c('L1','L2','tau3','tau4')],query = sample_LM,k = 10)$nn.idx
b = data.frame(Lcv = sample_LM[2]/sample_LM[1], tau3 = sample_LM[3], tau4 = sample_LM[4])
I = RANN::nn2(a[,c('Lcv','tau3','tau4')],query = b,k = 10)$nn.idx
# I = which.min(abs((data[,'L1']-sample_LM[,'L1'])/(sample_LM[,'L1']+0.0001)) +
#                 abs((data[,'L2']-sample_LM[,'L2'])/(sample_LM[,'L2']+0.0001)) +
#                 abs((data[,'tau3']-sample_LM[,'tau3'])/(sample_LM[,'tau3']+0.0001)) +
#                 abs((data[,'tau4']-sample_LM[,'tau4'])/(sample_LM[,'tau4']+0.0001)))
init_scale = Burr_InitValues[I,"init_scale"]
init_shape1 = Burr_InitValues[I,"init_shape1"]
init_shape2 = Burr_InitValues[I,"init_shape2"]
# init_scale = predict(scale_model, sample_LM)
# init_shape1 = predict(shape1_model, cbind(init_scale,sample_LM))
# init_shape2 = predict(shape2_model, cbind(init_scale,init_shape1,sample_LM))
start = c(init_scale, init_shape1, init_shape2)

para = pelp(lmom = lm[1:3],
            pfunc = pfunction,
            start = Burr_InitValues[I[1],c(1:3)],
            bounds = c(0, Inf),
            type = 's')
params = as.list(para$para)
params$location = loc

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
