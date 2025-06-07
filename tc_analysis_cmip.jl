#-----------------------------------------------------------------------------
# tc_analysis_cmip.jl
#
# read in SST (tos) and MPI(usually from scripts run on Casper)
# compute a timeseries of RONI from the SST
# compute the relative SST
#
# plot the roni ts, the relative SST, and the MPI
#
# to run:  
# >julia tc_analysis_cmip.jl
#
# we need to output the MPI composite data as a netcdf file so that it can be 
# averaged with the output from other models easily. 
#
# levi silvers                                            may 2025
#-----------------------------------------------------------------------------

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics
#using GLM
using DataFrames

include("ensoFuncs.jl")

lpath="/Users/C823281551/data/cmip6/"

#file1   = lpath*"MPIESM/potInt/tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012.73x144.nc"
#file1b  = lpath*"MPIESM/potInt/tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012_summer.73x144.nc"
#file2   = lpath*"MPIESM/potInt/MPI_ESM1_full_output.nc"
#file1   = lpath*"ACCESS-CM2/tos_73x144_pp.nc"
#file1b  = lpath*"ACCESS-CM2/tos_73x144_pp_summer.nc"
#file2   = lpath*"ACCESS-CM2/ACCESS_CM2_full_output.nc"
#file1   = lpath*"MIROC6/tos_73x144_pp.nc"
#file1b  = lpath*"MIROC6/tos_73x144_pp_summer.nc"
#file2   = lpath*"MIROC6/MIROC6_full_output.nc"
file1   = lpath*"MRI-ESM2-0/tos_73x144_pp.nc"
file1b  = lpath*"MRI-ESM2-0/tos_73x144_pp_summer.nc"
file2   = lpath*"MRI-ESM2-0/MPI_MRI-ESM2_full_output.nc"

#tag = "MPI_ESM"
#tag = "ACCESS_CM2"
#tag = "MIROC6"
tag = "MRI_ESM2_pylis"

inpFile = file1 

# define an array that mirrors the time period of the data
C = collect(2015.083333:1/12:2101)

# these will need to be adjusted for each, stupid, pathetic, model:
ensoDef   = 1.5; # MPIESM used 1.6 threshold, ACCESS-CM2 used 1.3
thresh = string(ensoDef)
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

sst_var    = data1["tos"]
sst_summer = data1b["tos"]
dims       = size(sst_var)
dims2      = size(sst_summer)
# grab the summer variables
#sst_var_s_index = Array{Union{Missing, Int}, 1}(undef, timelenH)
#sst_summer      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], timelenH)
#global ix = 1
##ix = 1
#for i in 5:12:timelen2
#    #println(i," ",i+1," ",i+2," ",i+3," ",i+4," ",i+5)
#    for j = 1:6
#        mnth = i+j
#        sst_var_s_index[ix]=mnth
#        global ix = ix + 1
#        #ix = ix + 1
#        #println("ix = ",ix)
#    end
#end
# pick out june-november from the full data set
# this line of code took 24 minutes to run...  unacceptable.
###sst_summer = sst_var[:,:,sst_var_s_index]
# as an alternative use the nco ncrcat function to subcycle through the full dataset:
#ncrcat -d time,6,,12,6 tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012.regrid.nc -O tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012_summer.nc
# use dimension time, starting at index 6, going to the end of the record, every 12 indices, keep the six 
# sequential indices.

sst_tr_mean          = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
timeAxis = collect(1.:1:timelenH); 

for i in 1:dims2[3]
    sst_tr_mean[i] = mean(skipmissing(sst_summer[:,tropS:tropN,i]))
end

agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
#rel_nina  = Array{Union{Missing, Float64}, 2}(undef, dims[1], 61)
#rel_nino  = Array{Union{Missing, Float64}, 2}(undef, dims[1], 61)
# loop over longitude and latitude:
# find the linear best fit regression at every grid-point
for i in 1:dims[1]
    for j in 1:dims[2]
    #for j in 60:120
        #agrid[i,j],bgrid[i,j] = find_best_fit(sst_tr_mean[:],sst_summer[i,j,:])
        agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis[:],sst_summer[i,j,:])
    end
end
#for i in 1:dims[1]
#    j0=1
#    for j in 60:120
#        rel_nina[i,j0],bgrid[i,j0] = find_best_fit(sst_tr_nina[:],sst_nina[i,j,:])  
#        rel_nino[i,j0],bgrid[i,j0] = find_best_fit(sst_tr_nino[:],sst_nino[i,j,:])
#        j0 = j0+1
#    end
#end

#timeAxis = collect(1.:1:timelenH); 
# timeAxis should be a monthly time series, but with only 6 months per year.
# find slope of tropical mean sst, in celcius per year
a,b = find_best_fit(timeAxis,sst_tr_mean)
slopePerCentury = a*6*100
relSST = agrid*6*100 .- slopePerCentury   

println("**************************************")
#print("max value is: ",maximum(relSST, dims = 2))
#print("min value is: ",minimum(relSST, dims = 2))
print("max value is: ",maximum(skipmissing(relSST)))
print("min value is: ",minimum(skipmissing(relSST)))
println("**************************************")
print("slot Per Century is: ",slopePerCentury)
println("**************************************")

# work on the MPI calculation
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
prepare_cmip_ts(inpFile,timelen2,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# use a running mean filter to smooth the timeseries.
smooth_12_ts(ts_roni,timelen2)
ts_roni_sm = ts_12_sm

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2];

# check thresholds for ENSO phases: 
thshd = ensoDef
check_thresh_high(inpFile, ts_roni_sm, thshd) # output is 'high[]'
thshd = -ensoDef
#thshd = -2.0
check_thresh_low(inpFile, ts_roni_sm, thshd) # output is 'low[]'

# create arrays for the composites of the PI
dims3     = size(vmax)
pi_high  = Array{Union{Missing, Float64}, 3}(undef, dims3[1], dims3[2], numfields)
pi_low   = Array{Union{Missing, Float64}, 3}(undef, dims3[1], dims3[2], numfields)

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,low[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,high[i]]   # El Nino months
end

blackbird = pi_high .- pi_low;
PI_comp_mn = mean(blackbird, dims=3)

#----------------------------------------------------------------------------------
## work on plots

levs = range(-10, 10, length = 21)
titsuf = ": Nino - Nina"
#tit = tag * titsuf
tit = "MPI:" * titsuf
fig = Figure(;
    size = (800,600),
    )
    ax = GeoAxis(fig[2,1];
      xticks = -180:30:180,
      yticks = -90:30:90,
      limits=(-180,180,-40,40),
      #xlabel = "time (yr)",
      #ylabel = "RONI",
      title  = "Relative SST trend (K/Century)")
      #fig_1_plot(PI_5_tmn[:,:,1],lon,lat,tit1,level1)
      bb = contourf!(ax, lon1d, lat1d, relSST[:,:],
               levels = range(-2, 2, length = 21), # rh
               colormap = :vik,
               extendlow = :auto, extendhigh = :auto
          )
    lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=1.25)
    Colorbar(fig[2,2], bb)
    ax = Axis(fig[1,1];
      xlabel = "time (yr)",
      ylabel = "RONI",
      #title  = tag * " threshold: " * ensoDef)
      title  = tag * " ENSO Threshold: " * thresh)
      #title  = "1 cmip6 model")
    lines!(C,ts_roni_sm, linestyle = :solid, linewidth=1.25)
    lines!(C,ensoPlotp, color = :black, linestyle = :solid, linewidth=1.25)
    lines!(C,ensoPlotn, color = :black,  linestyle = :solid, linewidth=1.25)
#plot!(C,ts_roni_sm) # plots the time series using dots for the data points. 
    plot!(C,forNew)
ax = GeoAxis(fig[3,1];
    xticks = -180:30:180,
    yticks = -90:30:90,
    limits=(-180,180,-40,40),
    xlabel = "time (yr)",
    ylabel = "RONI",
    title  = tit)
    bb = contourf!(ax, lon, lat, PI_comp_mn[:,:,1],
             levels = levs,
             colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
    lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=1.25)
    Colorbar(fig[3,2], bb)
fig


figname=tag*"_mpi.png"
save(figname, fig)


#----------------------------------------------------------------------------
## save output to a netcdf file
ds = NCDataset("test.nc","c")

defDim(ds,"lon",144)
defDim(ds,"lat",73)
defDim(ds,"time",1)

ds.attrib["title"] = "this is a test file"

#v1 = defVar(ds,"MPI",Float32,("lon","lat","time"))
v2 = defVar(ds,"lat1",Float32,("lat"))
#v3 = defVar(ds,"lon",Float32,("lon"))

PI_comp_mn_nan = nomissing(PI_comp_mn,NaN)

#v1[:,:,:] = PI_comp_mn_nan

v2[:]     = lat
#v3[:]     = lon

v1.attrib["units"] = "meters per second"
v1.attrib["commments"] = "this is the time avegage maximum potential intensity for TCs"
#v1.attrib["_FillValue"] = "Missing"

close(ds)
