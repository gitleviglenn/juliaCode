# plot linear sst trend from October and November over 1979-2024
#
# levi silvers

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics
using GLM
using Polynomials

include("ensoFuncs.jl")

path="/Users/C823281551/data/ERA5/"

filein1 = path*"era5_sst_1979th2024_OctNov_360x180.nc"
filein  = path*"MPI_ERA5_full_output.nc"
tag = "ERA5"
data1  = NCDataset(filein1)
data   = NCDataset(filein)

#lat = data["latitude"]
#lon = data["longitude"]
lat = data["lat"]
lon = data["lon"]
tme = data["valid_time"]
 
dim_var = data1["sst"]

timeAxis1 = collect(1:1:92);

#sst_var = data["sst"];
sst_var = data["vmax"];
dims    = size(dim_var)

# vmax has 408 timesteps.... sst_var has 92
# 10,11,22,23
#
# save i and i+1, jump 12, do again
for i in 10:12:408
  println("hah!")
  #j = i 
  #println(j)
  #j = j+1
  #println(j)
  #j = j+12
  println(i)
  i = i+1
  println(i)
end 

print("size of sst Array is: ",size(dim_var))
print("size of PI Array is: ",size(sst_var))

agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
# loop over longitude and latitude:
for i in 1:dims[1]
    for j in 1:dims[2]
        agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis1[:],sst_var[i,j,:])
        #agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis1[:],sst_var[:,i,j])
    end
end

levs = range(-.5, .5, length = 21)
blah = fig_anom_plot(agrid.*20,lon,lat,"linear trends October-November (C/decade)",levs)

save("era5_PI_LinTrend_OctThNov_1979th2024.png", blah, px_per_unit=6.0)
