#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot2Dearth.jl
#
# simple script to plot 2D data on a geohysical Earth grid
#
# https://geo.makie.org/v0.7.4/introduction
#
# plots continent lines
# changes the default tickmarks for the lat and lon grid that is shown
#
# colortables: 
#     batlow: a sequential gradient colortable that inclues beach, tan, green, 
#             and deep green/blue
#     vik:    a diverging gradient colortable that includes blues and reds (white)
#             in the middle. 
# 
# what are the junk noisy white areas in the GeoMakie projections?  
#
# ] add ColorSchemes
#
# https://juliadatascience.io/makie_colors
# https://docs.makie.org/dev/explanations/colors
# 
# the julia data science document I have been relying on doesn't go into 
# geophysical plotting
#
# levi silvers                                                oct 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics
using ColorSchemes
#using Plots

filein  = "/Users/C823281551/data/ERA5/era5_nina_2d.nc"
file2in = "/Users/C823281551/data/ERA5/era5mon1940toPresent_div_U_V.nc"

# read data from the incoming netcdf file: 
data  = NCDataset(filein)
data2 = NCDataset(file2in)

data.attrib

sst = data["sst"]
lat = data["latitude"]
lon = data["longitude"]
tim = data["time"]

velu = data2["u"]
velv = data2["v"]
div  = data2["d"]

# computing the average over all times is very slow.   why?
#sst_mn = mean(sst, dims=3)

sst_1t = sst[:,:,30]
u_1    = velu[:,:,1,500]
v_1    = velv[:,:,1,500]
d_1    = div[:,:,1,500]

  println("~~~~obs nino3.4~~~~~")
  println(collect(sst[600, :, 10]))
  println("~~~~obs trop mn~~~~~")

sst_tmn = mean(filter(!isnan, skipmissing(u_1)), dims=3)

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  isnan(sst[23,23,23])
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#if isnan(sst_1t)
#    println("not a number!")
#else
#    println("sst_1t is a number")
#end

#println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#println(lon[10:61])
#println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# for documentation on GeoMakie, see: 
# https://github.com/MakieOrg/GeoMakie.jl

sst_a = reshape(sst_tmn, (1440, 721))


function fig_1_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(600,300),
        )
    ax = Axis(f2[1,1];
        xticks = -180:30:180, 
        yticks = -90:30:90,
        xlabel="latitude",
        ylabel="longitude",
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             levels = range(-20, 20, length = 100), 
             #colormap = :batlow,
             colormap = :vik,
        )
        Colorbar(f2[1,2], bb)
    return f2
end

function fig_2_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(600,300),
        )
    ax = GeoAxis(f2[1,1]; 
        xticks = -180:30:180, 
        yticks = -90:30:90,
        xlabel="latitude",
        ylabel="longitude",
        title=tit,
        )
        bb   = contourf!(ax, d1, d2, inpv, 
               levels = range(-10, 10, length = 100),
               extendlow = :white, extendhigh = :white
        )
        Colorbar(f2[1,2], bb)
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
    return f2
end

# trying to shift longitude... still not fully working
#lon_0 = 90.0
#(blah_lon, blah_field) = cshift(lon, v_1, lon_0)

#fig1 = fig_1_plot(sst_a,lon,lat,"SST")
fig2 = fig_1_plot(v_1,lon,lat,"velocity") # uses equidistant.    
#fig2 = fig_2_plot(v_1,lon,lat,"velocity") # uses GeoAxis, and doesn't yet look good. 






