using NetCDF
using CairoMakie
# may need to install PyCall and PyPlot

# using Plots # --> should this be used as well?

# run with: 
# /Users/silvers/.juliaup/bin/julia era5map.jl

file1 = "/Users/silvers/data/ERA5/omega500_monthly_1979th2021_era5.nc"
file2 = "/Users/silvers/data/ERA5/era5_monthly_1959present_10mWindSpeed.nc"

# to print meta data to screen:
#ncinfo(file1)

# save the var w to test1:
test1=ncread(file1, "w")
test2=ncread(file2, "si10")

# how do I print out info about test1?  like max and min?
size(test1)
lons = ncread(file1, "longitude")
lats = ncread(file1, "latitude")

# plot(heatmap(x=lons,y=lats,z=test1))

fig = Figure()

# this figure looks reasonable: 
#fig = contour(lons,lats,test1[:,:,1])

#ax = Axis(fig[1,1],
Axis(fig[1,1],
  title = "in the aeroplane over the sea",
  xlabel = "up and over",
  ylabel = "no" 
)

#contour!(lons,lats,test1[:,:,1])
contour!(lons,lats,test2[:,:,12])

# this figure looks wrong:  why?
#contour(lons,lats,test1[:,:,30])

save("myfigureWind3.png", fig)

