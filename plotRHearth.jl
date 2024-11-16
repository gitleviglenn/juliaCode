#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot2Dearth.jl
#
# make a contour plot of mid-tropospheric RH
# 
# levi silvers                                                 nov 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics
using ColorSchemes
#using Plots

filein  = "/Users/C823281551/data/cmip6/CESM2/hur_Amon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215.nc"
#filein  = "/Users/C823281551/data/ERA5/era5_nina_2d.nc"
#file2in = "/Users/C823281551/data/ERA5/era5mon1940toPresent_div_U_V.nc"

# read data from the incoming netcdf file: 
data  = NCDataset(filein)
#data2 = NCDataset(file2in)

data.attrib

rh   = data["hur"] # hur(time, plev, lat, lon)
lat  = data["lat"]
lon  = data["lon"]
tme = data["time"]
lev  = data["plev"]

tim = 5
#rh1 = rh[:,tim,1,:]
#rh1 = rh[tim,1,:,:]
#rh1 = rh[:,:,tim,1]
rh1 = rh[:,:,1,tim]
#rh2 = rh[:,tim,2,:]
rh2 = rh[:,:,2,tim]

#velu = data2["u"]
#velv = data2["v"]
#div  = data2["d"]

function fig_1_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(600,300),
        )
    ax = Axis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             levels = range(0, 100, length = 10), 
             colormap = :Blues_8,
             #colormap = :navia,
             #colormap = :batlow,
             #colormap = :vik,
        )
        Colorbar(f2[1,2], bb)
    return f2
end
function fig_diff(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(600,300),
        )
    ax = Axis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             levels = range(-20, 20, length = 20), 
             #colormap = :Blues_8,
             #colormap = :navia,
             #colormap = :batlow,
             colormap = :vik,
        )
        Colorbar(f2[1,2], bb)
    return f2
end

fig = fig_1_plot(rh2,lon,lat,"joe chr")



