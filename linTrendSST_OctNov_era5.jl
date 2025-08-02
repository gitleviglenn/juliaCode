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

filein  = path*"era5_sst_1979th2024_OctNov_360x180.nc"
tag = "ERA5"
data   = NCDataset(filein)

#lat = data["latitude"]
#lon = data["longitude"]
lat = data["lat"]
lon = data["lon"]
tme = data["valid_time"]
 
sst_var = data["sst"]

timeAxis1 = collect(1:1:92);

sst_var = data["sst"];
dims    = size(sst_var)

agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
# loop over longitude and latitude:
for i in 1:dims[1]
    for j in 1:dims[2]
        agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis1[:],sst_var[i,j,:])
    end
end

levs = range(-.5, .5, length = 21)
blah = fig_anom_plot(agrid.*20,lon,lat,"linear trends October-November (C/decade)",levs)

save("era5_sst_LinTrend_OctThNov_1979th2024.png", blah, px_per_unit=6.0)
