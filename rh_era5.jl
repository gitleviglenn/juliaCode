#-----------------------------------------------------------------------------------------------
# rh_era5.jl
#
#-----------------------------------------------------------------------------------------------

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics

include("ensoFuncs.jl")
#
tag = "Dude"
# create the indices that correspond to Nino and Nina years/months
#fig2name = tag*"_rh_SH_test.png"
# NH
fig2name = tag*"_rh_nino_comp_NH_test.png"
ninoyears = [18 54 90 150 174 234 306 402]
ninayears = [102 114 210 246 318 366 378 390]
# SH
#fig2name = tag*"_rh_nino_comp_SH_test.png"
#ninoyears = [23 35 59 96 155 239 311 347]
#ninayears = [107 119 215 251 263 335 371 383]
function create_indices(years)
  ensoInd = Matrix{Int64}(undef, 8, 6)
  for i in 1:8
    si = years[i]
    ensoInd[i,:] = si:si+5
  end
  #o1=years[1]:years[1]+5
  #ensoInd = [o1 o2 o3 o4 o5 o6 o7 o8]
  return ensoInd
end 

ninoInd = create_indices(ninoyears)
ninaInd = create_indices(ninayears)
#

path="/Users/C823281551/data/ERA5/"

filein  = path*"era5_rh_1990th2023_360x80.nc"
tag = "ERA5"
data   = NCDataset(filein)

lat = data["lat"]
lon = data["lon"]
lev = data["pressure_level"]

rh_var = data["r"]

numfields = 48 
numall    = 408

lat1      = 1
lat2      = 81

dims = size(rh_var)

rh_1            = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
rh_2            = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
rh_low          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
rh_high         = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
rh_tot          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numall)

endi = 48
# low is an array that contains the timesteps representing the negative phase of ENSO
# high is an array that contains the timesteps representing the positive phase of ENSO
low = ninaInd
high = ninoInd

# select the case determining if rh will be plotted at one level or at the average of two
pcase = 0# pcase = 0 corresponds to plotting rh on 1 level
# calculate total rh field, for all times
level=2 # level 2 should correspond to the 700 hPa pressure level. 
if pcase == 0
  for i in 1:408
    rh_tot[:,lat1:lat2,i]   = rh_var[:,lat1:lat2,level,i] 
  end
  for i in 1:endi
    # low values, la nina
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,level,low[i]]
    rh_low[:,lat1:lat2,i]   = rh_2[:,lat1:lat2,1,i] 
    # high values, el nino
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,level,high[i]]
    rh_high[:,lat1:lat2,i]  = rh_2[:,lat1:lat2,1,i] 
  end
else
  for i in 1:408
    rh_tot[:,lat1:lat2,i]   = (rh_var[:,lat1:lat2,2,i] + rh_var[:,lat1:lat2,1,i])./2
  end
  for i in 1:endi
    # low values, la nina
    rh_1[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,1,low[i]]
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,2,low[i]]
    rh_low[:,lat1:lat2,i]   = (rh_2[:,lat1:lat2,1,i] + rh_1[:,lat1:lat2,1,i])./2
    # high values, el nino
    rh_1[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,1,high[i]]
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,2,high[i]]
    rh_high[:,lat1:lat2,i]  = (rh_2[:,lat1:lat2,1,i] + rh_1[:,lat1:lat2,1,i])./2
  end
end 


rh_tot_tmn  = mean(rh_tot, dims=3)
rh_high_tmn = mean(rh_high, dims=3)
rh_low_tmn  = mean(rh_low, dims=3)

data_2_plot_tot  = rh_tot_tmn
data_2_plot_anom = rh_high_tmn - rh_low_tmn

tit="ERA5 RH Composite (%)"
levs = range(-15., 15., length = 21)
#fig2name = tag*"_rh_nino_comp_SHmn.png"
#fig = fig_tot_plot(data_2_plot_tot[:,:,1],lon,lat,tit)
fig = fig_anom_plot(data_2_plot_anom[:,:,1],lon,lat,tit,levs)
#lonlon=[1:360;]

#fig = fig_anom_plot(data_2_plot_anom[:,:],lon,lat,tit)
save(fig2name, fig)
