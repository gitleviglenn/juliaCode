#-----------------------------------------------------------------------------------------------
# vShearMap.jl
#-----------------------------------------------------------------------------------------------
using CairoMakie
using GeoMakie
using NCDatasets

path="/Users/C823281551/data/"
modelp="MPI-ESM1-2-LR" 
#filehur  = path*"cmip6/MPIESM/hur_Amon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216.nc" 

#filein  = path*"cmip6/MPIESM/hur_Amon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_regridded3.nc" 
#filein  = path*"cmip6/MPIESM/tos_Omon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_pm40b.nc" 

filein   = path*"cmip6/CESM2/ua_Amon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215_360x180.nc"
filein2  = path*"cmip6/CESM2/va_Amon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215_360x180.nc"

tag = "CESM2"

data   = NCDataset(filein)
data2  = NCDataset(filein2)

data.attrib

lat  = data["lat"]
lon  = data["lon"]
lev  = data["plev"]
tme = data["time"]

u_var  = data["ua"] # hur(time, plev, lat, lon)
u_var  = data["ua"] # hur(time, plev, lat, lon)
v_var  = data2["va"] # hur(time, plev, lat, lon)
v_var  = data2["va"] # hur(time, plev, lat, lon)
tit = "wind: CESM2"

# tropics in CESM2
#lat1=50
#lat2=130
lat1=1
lat2=180
plev=2

tim = 506 # time step at which data will be plotted.

data_hi  = u_var[:,lat1:lat2,2,tim]
data_low = u_var[:,lat1:lat2,1,tim]
uSh      = data_hi - data_low
data_hi  = v_var[:,lat1:lat2,2,tim]
data_low = v_var[:,lat1:lat2,1,tim]
vSh      = data_hi - data_low

dims = size(u_var)
u_2  = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
u_1  = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
uSh_tmp  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
v_2  = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
v_1  = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
vSh_tmp  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
VWS_low_full  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
VWS_high_full  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
u_2b  = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
u_1b  = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
uSh_tmpb  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
v_2b  = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
v_1b  = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
vSh_tmpb  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)

endi = numfields
# for low values of ENSO timeseries.   Also compute for high values.  
for i in 1:endi
  u_2[:,lat1:lat2,1,i] = u_var[:,lat1:lat2,2,low[i]]
  u_1[:,lat1:lat2,1,i] = u_var[:,lat1:lat2,1,low[i]]
  uSh_tmp[:,lat1:lat2,i] = u_2[:,lat1:lat2,1,i] - u_1[:,lat1:lat2,1,i]
  v_2[:,lat1:lat2,1,i] = v_var[:,lat1:lat2,2,low[i]]
  v_1[:,lat1:lat2,1,i] = v_var[:,lat1:lat2,1,low[i]]
  vSh_tmp[:,lat1:lat2,i] = v_2[:,lat1:lat2,1,i] - v_1[:,lat1:lat2,1,i]
  VWS_low_full[:,lat1:lat2,i] = sqrt.(uSh_tmp[:,lat1:lat2,i].^2 .+ vSh_tmp[:,lat1:lat2,i].^2)
  #vSh[:,lat1:lat2,i] = v_2 - v_1
  u_2b[:,lat1:lat2,1,i] = u_var[:,lat1:lat2,2,high[i]]
  u_1b[:,lat1:lat2,1,i] = u_var[:,lat1:lat2,1,high[i]]
  uSh_tmpb[:,lat1:lat2,i] = u_2b[:,lat1:lat2,1,i] - u_1b[:,lat1:lat2,1,i]
  v_2b[:,lat1:lat2,1,i] = v_var[:,lat1:lat2,2,high[i]]
  v_1b[:,lat1:lat2,1,i] = v_var[:,lat1:lat2,1,high[i]]
  vSh_tmpb[:,lat1:lat2,i] = v_2b[:,lat1:lat2,1,i] - v_1b[:,lat1:lat2,1,i]
  VWS_high_full[:,lat1:lat2,i] = sqrt.(uSh_tmpb[:,lat1:lat2,i].^2 .+ vSh_tmpb[:,lat1:lat2,i].^2)
end

# need to compute the average of uSh_tmp and vSh_tmp over all i values
# do the same for the high values, and then take the difference.  

VWS_low  = mean(VWS_low_full, dims = 3)
#vSh_low  = mean(vSh_low_full, dims = 3)
VWS_high = mean(VWS_high_full, dims = 3)
#vSh_high = mean(vSh_high_full, dims = 3)

# should total shear be computed before averagine over time? 
#uSh = 

data_2_plot = VWS_high - VWS_low
#data_2_plot = sqrt.(uSh.^2 .+ vSh.^2)
#data_2_plot = sqrt.(uSh.^2)

function fig_1_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(600,300),
        )
    #ax = Axis(f2[1,1]; #--> default plot is rectangular equidistant 
    ax = GeoAxis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        limits=(-180,180,-40,40),
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             #levels = range(0, 50, length = 25), # tos
             levels = range(-15, 15, length = 20), # rh
             #colormap = :Blues_8,
             #colormap = :navia,
             #colormap = :batlow,
             colormap = :vik,
        )
        lines!(ax, GeoMakie.coastlines())
        Colorbar(f2[1,2], bb)
    return f2
end

fig = fig_1_plot(data_2_plot[:,:,1],lon,lat[lat1:lat2],tit)
filename="vertShear"*tag*"comp.png"
save(filename, fig)
#save("vertShear.png", fig)

