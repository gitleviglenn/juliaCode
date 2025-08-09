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

function local_fig(inpv,d1,d2,tit,levs)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        #figure_padding=(10,15,10,10),
        backgroundcolor=:white,
        size=(600,400),
        #size=(900,400),
        #size=(600,300), # this increases tickfont size, but doesn't print -30S!!#$%
        )   
    ax = GeoAxis(f2[1,1];
        dest="+proj=latlon",
        #xticks = -180:30:180,
        xticks = -110:10:-10, 
        #yticks = -90:30:90,
        yticks = 0:10:30,
        xlabel="longitude",
        ylabel="latitude",
        #limits=(-180,180,-90,90),
        limits=(-110,-10,0,30),
        title=tit,
        #xticklabelsize = 22, # 14,16 are pretty reasonable sizes
        #yticklabelsize = 22, # 22 used for 8 panel figure that needs larger font
        xticklabelsize = 14, # 14,16 are pretty reasonable sizes
        yticklabelsize = 14, # 22 used for 8 panel figure that needs larger font
        )   
        bb = contourf!(ax, d1, d2, inpv,
             levels = levs,
             #colormap = :batlow,
             #colormap = :bam, # default for shear plot (greens and pinks)
             #colormap = :seismic, # colors are a bit harsh
             colormap = :vik, # default for redish bluish for relative SST
             #colormap = :BrBg, # better for RH  browns and greens
             #colormap = :roma,
             extendlow = :auto, extendhigh = :auto
        )   
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[1,2], bb) 
    return f2
end

path="/Users/C823281551/data/ERA5/"

filein1  = path*"era5_sst_1979th2024_OctNov_360x180.nc"
filein   = path*"MPI_ERA5_OctNov_full_output.nc"
filein2  = path*"era5_hur_OctNov_1979th2024_360x180.nc"
filein3  = path*"era5_uWind_OctNov_1979th2024_360x180.nc"
filein4  = path*"era5_vWind_OctNov_1979th2024_360x180.nc"

tag = "ERA5"
data1  = NCDataset(filein1)
data2  = NCDataset(filein2)
data3  = NCDataset(filein3)
data4  = NCDataset(filein4)
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
rh_var  = data2["r"];
u_var   = data3["u"];
v_var   = data4["v"];

# relative humidity
level=1 # level 2 should correspond to the 700 hPa pressure level. 
lat1 = 1
lat2 = 180

rh_lev       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], dims[3])
u_tot        = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], dims[3])
v_tot        = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], dims[3])
VWS_tot      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], dims[3])

println("size of rh_lev is: ",size(rh_lev))

for i in 1:dims[3]
  #rsst_test      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 1)
  rh_lev[:,lat1:lat2,i]   = rh_var[:,lat1:lat2,level,i]
  #rh_lev[:,:,i] = rh_test # ?? should be the average of 6 fields...mpi_b[:,:,]
  u_tot[:,lat1:lat2,i]  = u_var[:,lat1:lat2,2,i] - u_var[:,lat1:lat2,1,i]
  v_tot[:,lat1:lat2,i]  = v_var[:,lat1:lat2,2,i] - v_var[:,lat1:lat2,1,i]
  VWS_tot[:,lat1:lat2,i]  = sqrt.(u_tot[:,lat1:lat2,i].^2 .+ v_tot[:,lat1:lat2,i].^2)
end

#  octNovInd = Matrix{Int32}(undef, 92)
octNovInd = Array{Union{Missing, Float64}, 1}(undef, dims[3])
pi_var = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], dims[3])
#
## save i and i+1, jump 12, do again
#for i in 10:12:408
#  if i==10 
#    global j=1
#    println("nope!")
#  end
#  octNovInd[j] = i
#  println(j)
#  j = j+1
#  octNovInd[j] = i+1
#  println(j)
#  j = j+1
#end 
#
#print(octNovInd[:])
#pi_var = sst_var[:,:,octNovInd]

print("size of sst Array is: ",size(dim_var))
print("size of PI Array is: ",size(sst_var))

#print(octNovInd[:])

agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
# loop over longitude and latitude:
for i in 1:dims[1]
    for j in 1:dims[2]
        #agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis1[:],pi_var[i,j,:])
        agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis1[:],sst_var[i,j,:])
    end
end

#levs = range(-.5, .5, length = 21)
levs = range(-2.5, 2.5, length = 21)
blah = local_fig(agrid.*20,lon,lat,"PI linear trends October-November (m/s)/decade",levs)
save("era5_PI_LinTrend_OctThNov_1979th2024_region.png", blah, px_per_unit=6.0)

agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
for i in 1:dims[1]
    for j in 1:dims[2]
        agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis1[:],rh_lev[i,j,:])
    end
end
levs = range(-1.0, 1.0, length = 21)
blah = local_fig(agrid.*20,lon,lat,"RH linear trends October-November %/decade",levs)
save("era5_RH850_LinTrend_OctThNov_1979th2024_region.png", blah, px_per_unit=6.0)

agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
for i in 1:dims[1]
    for j in 1:dims[2]
        agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis1[:],VWS_tot[i,j,:])
    end
end
levs = range(-2.5, 2.5, length = 21)
blah = local_fig(agrid.*20,lon,lat,"VWS linear trends October-November (m/s)/decade",levs)
save("era5_VWS_LinTrend_OctThNov_1979th2024_region.png", blah, px_per_unit=6.0)







