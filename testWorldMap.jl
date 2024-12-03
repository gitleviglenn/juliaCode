#----------------------------------------------------------
# testWorldMap.jl
#
# to change the variable being plotted, one must change: 
#   filein
#   make sure the dimensions of data_2_plot are correct 
#   the variable details in the variable details block 
#   the 'levels' defined within fig_1_plot
#
# regridding the SST data can be tricky.   for CESM2 I used:
# cdo -L remapbil,mygrid tos_Omon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215.nc tos_Omon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215_regrid2.nc
#
# for MPI I used: 
# ncatted -a coordinates,tos,c,c,latitude longitude tos_Omon_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_20150116-21001216_pm40b.nc
# cdo -L remapbil,mygrid2 -sethalo,-1,-1 tos_Omon_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_20150116-21001216.nc tos_Omon_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_20150116-21001216_pm40b.nc
# according to Uwe's comment here: https://code.mpimet.mpg.de/boards/2/topics/15918?r=15935#message-15935
# data on a MPIOM model grid have 2 extra columns, which can be removed with the sethalo command
#
# levi silvers                                    nov 2024
#----------------------------------------------------------
using CairoMakie
using GeoMakie
using NCDatasets

path="/Users/C823281551/data/"
modelp="MPI-ESM1-2-LR" 
#filehur  = path*"cmip6/MPIESM/hur_Amon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216.nc" 

#filein  = path*"cmip6/MPIESM/hur_Amon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_regridded3.nc" 
#filein  = path*"cmip6/MPIESM/tos_Omon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_pm40b.nc" 

filein  = path*"cmip6/CESM2/tos_Omon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215_regrid2.nc"

data  = NCDataset(filein)

data.attrib

#var   = data["hur"] # hur(time, plev, lat, lon)
#var  = data["tos"] # hur(time, plev, lat, lon)
lat  = data["lat"]
lon  = data["lon"]
tme = data["time"]
#lev  = data["plev"]

#----------------------------------------------------------
# variable details: 

## tos
var  = data["tos"] # hur(time, plev, lat, lon)
tit = "tos: CESM2 b"
# MPI
#lat1=1
#lat2=80
# CESM2
lat1=50
lat2=130

## rh
#var   = data["hur"] # hur(time, plev, lat, lon)
#tit = "rh"
#lat1 = 51
#lat2 = 130
#----------------------------------------------------------

tim = 506 # time step at which data will be plotted.
data_2_plot = var[:,lat1:lat2,tim]
#data_2_plot = var[:,lat1:lat2,2,tim]

#rh1 = rh[:,:,1,tim]
#rh2 = rh[:,:,2,tim]
# rh
#data_2_plot = var[:,51:130,2,tim]
# tos
#data_2_plot = var[:,:,tim]
#s = circshift(lon, 180) # passing this in as the lon dimension shifts the rh but not coastlines. 

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
        limits=(-180,180,-40,40),
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             levels = range(0, 40, length = 10), # tos
             #levels = range(0, 100, length = 10), # rh
             colormap = :Blues_8,
             #colormap = :navia,
             #colormap = :batlow,
             #colormap = :vik,
        )
        lines!(ax, GeoMakie.coastlines())
        Colorbar(f2[1,2], bb)
    return f2
end

# rh
#data_2_plot = var[:,51:130,2,tim]
# tos
#data_2_plot = var[:,:,tim]
#s = circshift(lon, 180) # passing this in as the lon dimension shifts the rh but not coastlines. 

#fig = fig_1_plot(data_2_plot,lon,lat[51:130],"relative humidity")
fig = fig_1_plot(data_2_plot,lon,lat[lat1:lat2],tit)
save("testMap.png", fig)

