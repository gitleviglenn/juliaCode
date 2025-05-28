# sst_mpi_cmip.jl
#
# no. 13
#
# produces a three panel plot showing the relative SST trend, the 
# relative MPI trend, and the climatological MPI.
#
# relative trends are computed with linear regression.   two regression 
# options are available, chosen with the varcase variable
#
# these panels were inspired by Vecchi and Soden, 2007, Effect of remote 
# sea surface temperature change on tropical cyclone potential intensity, nature
#
# levi silvers                              may 26th, 2025

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics
#using GLM
using DataFrames

include("ensoFuncs.jl")

lpath="/Users/C823281551/data/cmip6/"

file1   = lpath*"MPIESM/potInt/tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012.73x144.nc"
file1b  = lpath*"MPIESM/potInt/tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012_summer.73x144.nc"
file2   = lpath*"MPIESM/potInt/MPI_ESM1_full_output.nc"
file2b  = lpath*"MPIESM/potInt/MPI_ESM1_full_summer.nc"
inpFile = file1 

# in some cases we will want to select the summer months from the MPI files.  for this we should be 
# able to use this nco command: 
#ncrcat -d time,6,,12,6 test.nc -O MPI_ESM1_full_summer.nc
# however, if time isn't the record variable, then this doesn't work.   it seems like to be the 
# record variable, time has to be unlimitted.   I don't understand this, but the ncks command
# below seems to fix the problem: 
#ncks --mk_rec_dmn time MPI_ESM1_full_output.nc -o test.nc

tag = "MPI_ESM"

# define an array that mirrors the time period of the data
C = collect(2015.083333:1/12:2101)

# these will need to be adjusted for each, stupid, pathetic, model:
ensoDef   = 1.6;
numfields = 60;
# these values are for a 360x180 grid
#lat1 = 71;lat2 = 110;lon34a = 10;lon34b = 61;lat34a = 85;lat34b = 96
# these values are for a 144x73 grid
lat1 = 29;lat2 = 45;lon34a = 5;lon34b = 25;lat34a = 35;lat34b = 39

timelen2 = 1032
timelenH = 516

# load the data
data1 = NCDataset(file1)
data1b= NCDataset(file1b)
data2 = NCDataset(file2)
data2b= NCDataset(file2b)
vmax  = data2["vmax"]
lat   = data2["lat"]
lon   = data2["lon"]
lat1d = data1b["lat"]
lon1d = data1b["lon"]

# these values only work on a 360x180 grid
#tropS = 50
#tropN = 130
# these values only work on a 144x73 grid
tropS = 20
tropN = 52

sst_var     = data1["tos"]
sst_summer  = data1b["tos"]
vmax_summer = data2b["vmax"]
dims        = size(sst_var)
dims2       = size(sst_summer)

sst_tr_mean          = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
timeAxis = collect(1.:1:timelenH); 

for i in 1:dims2[3]
    sst_tr_mean[i] = mean(skipmissing(sst_summer[:,tropS:tropN,i]))
end

agrid   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
agrid2  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid2   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
agridb  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgridb  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
agridb2 = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgridb2  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
#rel_nina  = Array{Union{Missing, Float64}, 2}(undef, dims[1], 61)
#rel_nino  = Array{Union{Missing, Float64}, 2}(undef, dims[1], 61)
# loop over longitude and latitude:
# find the linear best fit regression at every grid-point
for i in 1:dims[1]
    for j in 1:dims[2]
    #for j in 60:120
        agrid[i,j],bgrid[i,j]     = find_best_fit(timeAxis[:],sst_summer[i,j,:])
        agrid2[i,j],bgrid2[i,j]   = find_best_fit(sst_tr_mean[:],sst_summer[i,j,:])
        agridb[i,j],bgridb[i,j]   = find_best_fit(timeAxis[:],vmax_summer[i,j,:])
        agridb2[i,j],bgridb2[i,j] = find_best_fit(sst_tr_mean[:],vmax_summer[i,j,:])
    end
end

#timeAxis = collect(1.:1:timelenH); 
# timeAxis should be a monthly time series, but with only 6 months per year.
# find slope of tropical mean sst, in celcius per year
a,b = find_best_fit(timeAxis,sst_tr_mean)
slopePerCentury = a*6*100

relSST = agrid*6*100 .- slopePerCentury   
relMPI = agridb*6*100 .- slopePerCentury   

println("**************************************")
print("max value is: ",maximum(skipmissing(relSST)))
print("min value is: ",minimum(skipmissing(relSST)))
println("**************************************")
print("slot Per Century is: ",slopePerCentury)
println("**************************************")


PI_mn = mean(vmax, dims=3)
sst_summ_mn = mean(sst_summer, dims=3)

#----------------------------------------------------------------------------------
## work on plots

# decide what to plot, either the change per degree of tropical warming, or the change
# per century.   these cases are equivalent to figure 1 and figure 2 of 
# Vecchi and Soden, 2007

varcase = 1

if varcase > 1 
  tit1 = "Relative SST trend (K/Century)"
  field2plotA = relSST
  levs1 = range(-2, 2, length = 21)
  tit2 = "Relative MPI trend ((m/s)/Century)"
  field2plotB = relMPI
  levs2 = range(-20, 20, length = 21)
else
  tit1 = "SST change per degree tropical SST warming (K/K)"
  field2plotA = agrid2 
  levs1 = range(-2, 2, length = 21)
  #field2plotA = agrid2 .- 1 # normalize by a 1K warming
  #levs1 = range(-1, 1, length = 21)
  tit2 = "MPI change per degree tropical SST warming ((m/s)/K)"
  field2plotB = agridb2
  levs2 = range(-10, 10, length = 21)
end

levs = range(-10, 10, length = 21)
titsuf = "Maximum Potential Intensity (m/s)"
#tit = tag * titsuf
tit = titsuf
fig = Figure(;
    size = (800,600),
    )
    ax = GeoAxis(fig[1,1];
      xticks = -180:30:180,
      yticks = -90:30:90,
      limits=(-180,180,-40,40),
      title  = tit1)
      #title  = "Relative SST trend (K/Century)")
      bb = contourf!(ax, lon1d, lat1d, field2plotA[:,:],
               levels = levs1, 
               colormap = :vik,
               extendlow = :auto, extendhigh = :auto
          )
    lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=1.25)
    Colorbar(fig[1,2], bb)
ax = GeoAxis(fig[2,1];
    xticks = -180:30:180,
    yticks = -90:30:90,
    limits=(-180,180,-40,40),
      title  = tit2)
    bb = contourf!(ax, lon, lat, field2plotB[:,:],
    #bb = contourf!(ax, lon, lat, sst_summ_mn[:,:,1],
             levels = levs2, 
             #levels = range(0, 32, length = 33),
             colormap = :vik,
             #colormap = :batlow,
             extendlow = :auto, extendhigh = :auto
        )
    lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=1.25)
    Colorbar(fig[2,2], bb)
ax = GeoAxis(fig[3,1];
    xticks = -180:30:180,
    yticks = -90:30:90,
    limits=(-180,180,-40,40),
    xlabel = "time (yr)",
    ylabel = "MPI",
    title  = tit)
    bb = contourf!(ax, lon, lat, PI_mn[:,:,1],
             #levels = levs,
             levels = range(40, 120, length = 21), # rh
             colormap = :batlow,
             extendlow = :white, extendhigh = :auto
             #extendlow = :auto, extendhigh = :auto
        )
    lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=1.25)
    Colorbar(fig[3,2], bb)
fig


#figname=tag*"_sst_mpi_relSST.png"
figname=tag*"_sst_mpi_mean.png"
save(figname, fig)


