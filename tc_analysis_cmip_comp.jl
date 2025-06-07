#-----------------------------------------------------------------------------
# tc_analysis_cmip_temp.jl
#
# read in SST (tos) and MPI(usually from scripts run on Casper)
# compute a timeseries of RONI from the SST
# compute the relative SST
#
# plot the roni ts, the relative SST, and the MPI
#
# changes from tc_analysis_cmip.jl:  we now use the sst from within the 
# model_full_output_##.nc file, rather than a separate sst file.  as a result
# the input file to the function that calculates roni is not the full_output
# files, and that function now expects the sst variable name instead of the 
# tos variable name.  
# 
# to run:  
# >julia tc_analysis_cmip.jl
#
# The work of Vecchi and Soden, 2007 focused on summer time, or NH TC season (June-November)
# months.   We can do that.   But if we want to look at ENSO based composites it makes more
# sense to start with the complete time series and then pick out the positive and negative
# phases.  
#
# we need to output the MPI composite data as a netcdf file so that it can be 
# averaged with the output from other models easily. 
#
# the nco ncrcat function can be used to subcycle through the full dataset:
# ncrcat -d time,6,,12,6 tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012.regrid.nc -O tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012_summer.nc
# use dimension time, starting at index 6, going to the end of the record, every 12 indices, keep the six 
# sequential indices.
# ncks --mk_rec_dmn time filein.nc -O fileout.nc
# ncrcat -d time,6,,12,6 filein.nc -O fileout.nc
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

#-------------------------------------------------------------------------------------------------------------------
# specify file paths and names:

lpath="/Users/C823281551/data/cmip6/potInt/"


#file1   = lpath*"MPIESM/potInt/tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012.73x144.nc"
#file1b  = lpath*"MPIESM/potInt/tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012_summer.73x144.nc"
#file2   = lpath*"MPIESM/potInt/MPI_ESM1_full_output.nc"

##file1   = lpath*"ACCESS-CM2/tos_73x144_pp.nc"
##file1b  = lpath*"ACCESS-CM2/tos_73x144_pp_summer.nc"
#file1   = lpath*"ACCESS_CM2_full_output_summer.nc"
#file2   = lpath*"ACCESS_CM2_full_output_record.nc"

##file1   = lpath*"MIROC6/tos_73x144_pp.nc"
##file1b  = lpath*"MIROC6/tos_73x144_pp_summer.nc"
#file1   = lpath*"MIROC6_full_output_summer.nc"
#file2   = lpath*"MIROC6_full_output_record.nc"

# MRI-ESM2
#file1   = lpath*"MPI_MRI-ESM2_full_output_summer.nc"
#file2   = lpath*"MPI_MRI-ESM2_full_output_record.nc"

##file1   = lpath*"MRI-ESM2-0/tos_73x144_pp.nc"
##file1b  = lpath*"MRI-ESM2-0/tos_73x144_pp_summer.nc"

# MPI-ESM1-2
file1   = lpath*"MPI_ESM1_2_full_output_summer.nc"
file2   = lpath*"MPI_ESM1_2_full_output_record.nc"

##file2b  = "/Users/C823281551/data/MPI-ESM1/tos_Omon_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_20150116-21001216_latlon.nc"

# GFDL
#file1   = lpath*"GFDL-ESM4_full_output_nohalo_summer.nc"
#file2   = lpath*"GFDL-ESM4_full_output_nohalo_record.nc"

# IPSL-CM6A
#file1   = lpath*"MPI_IPSL-CM6A_full_output_summer.nc"
#file2   = lpath*"MPI_IPSL-CM6A_full_output_record.nc"

#-------------------------------------------------------------------------------------------------------------------
# determine which file will be used to compute RONI:
inpFile = file2 
#inpFile = file2b 

#-------------------------------------------------------------------------------------------------------------------
# tag to be used in plot title and figure name:

tag = "MPI_ESM1_2"
#tag = "ACCESS_CM2"
#tag = "MIROC6"
#tag = "MRI_ESM2"
#tag = "GFDL_ESM4"
#tag = "MPI_IPSL"

# define parameters
#
# these will need to be adjusted for each, stupid, pathetic, model:
ensoDef   = 1.0; # MPIESM used 1.6 threshold, ACCESS-CM2 used 1.3
#ensoDef   = 0.5; # MPIESM used 1.6 threshold, ACCESS-CM2 used 1.3
thresh = string(ensoDef)
numOfevents = 10;
numfields = 6*numOfevents;
# these values are for a 360x180 grid
#lat1 = 71;lat2 = 110;lon34a = 10;lon34b = 61;lat34a = 85;lat34b = 96
# these values are for a 144x73 grid
lat1 = 29;lat2 = 45;lon34a = 5;lon34b = 25;lat34a = 35;lat34b = 39
# below are indices for a 360x40 degree resolution file from MPI-ESM: testing only
#lat1 = 1;lat2 = 40;lon34a = 10;lon34b = 61;lat34a = 15;lat34b = 25

timelen      = 1032
timelenH     = 516
timeAxis     = collect(1.:1:timelenH); 
timeAxisFull = collect(1.:1:timelen); 
# define an array that mirrors the time period of the data
C = collect(2015.083333:1/12:2101)

# load the data
#data1 = NCDataset(file1)
data1 = NCDataset(file1)
data2 = NCDataset(file2)
#data2b = NCDataset(file2b)
vmax  = data2["vmax"]
lat   = data2["lat"]
lon   = data2["lon"]
lat1d = data1["lat"]
lon1d = data1["lon"]

# these values only work on a 360x180 grid
#tropS = 50
#tropN = 130
# these values only work on a 144x73 grid
tropS = 21
tropN = 53

println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1")
println(lat1d[tropS:tropN])
println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1")

sst_var     = data2["sst"]
sst_summer  = data1["sst"]
dims        = size(sst_var)
dims2       = size(sst_summer)

#dims = size(vmax)

sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
#rel_nina  = Array{Union{Missing, Float64}, 2}(undef, dims[1], 61)
#rel_nino  = Array{Union{Missing, Float64}, 2}(undef, dims[1], 61)

##-------------------------------------------------------------------------------
## work on the ralative change of SST using linear regression
#for i in 1:dims2[3]
#    sst_tr_mean[i] = mean(skipmissing(sst_summer[:,tropS:tropN,i]))
#end
#
## loop over longitude and latitude:
## find the linear best fit regression at every grid-point
#for i in 1:dims[1]
#    for j in 1:dims[2]
#    #for j in 60:120
#        #agrid[i,j],bgrid[i,j] = find_best_fit(sst_tr_mean[:],sst_summer[i,j,:])
#        agrid[i,j],bgrid[i,j]   = find_best_fit(timeAxis[:],sst_summer[i,j,:])
#        a_full[i,j],b_full[i,j] = find_best_fit(timeAxisFull[:],sst_var[i,j,:])
#    end
#end
#
## timeAxis should be a monthly time series, but with only 6 months per year.
## find slope of tropical mean sst, in celcius per year
#a,b = find_best_fit(timeAxis,sst_tr_mean)
#slopePerCentury = a*6*100
#relSST = agrid*6*100 .- slopePerCentury   
#
#println("**************************************")
#print("max value is: ",maximum(skipmissing(relSST)))
#print("min value is: ",minimum(skipmissing(relSST)))
#println("**************************************")
#print("slope Per Century is: ",slopePerCentury)
#println("**************************************")

#-------------------------------------------------------------------------------
# work on the MPI calculation
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn


# use a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm

# choose not to use the 12 point running mean: 
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slops that we want 
# to save for each model

# check thresholds for ENSO phases: 
## This method grabs a given number of point that are above the threshold
#thshd = ensoDef
#check_thresh_high(inpFile, ts_roni_sm, thshd) # output is 'high[]'
#thshd = -ensoDef
#check_thresh_low(inpFile, ts_roni_sm, thshd) # output is 'low[]'

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_indices_Nino=padding(local_max_indices)
pad_indices_Nina=padding(local_min_indices)

println("pad_indices_Nina: ",pad_indices_Nina[:])
println("size of vmax is: ",size(vmax))

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_indices_Nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_indices_Nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#-------------------------------------------------------------------------------
# work on the ralative change of SST using linear regression
for i in 1:dims2[3]
    sst_tr_mean[i] = mean(skipmissing(sst_summer[:,tropS:tropN,i]))
end

for i in 1:numfields
  enso_low[:,:,i]  = sst_var[:,:,pad_indices_Nina[i]]    # La Nina months
  enso_high[:,:,i] = sst_var[:,:,pad_indices_Nino[i]]    # El Nino months
end   

# I think that I need to compute the relative sst change for Nino and Nina 
# separately, then subtract the slopePerCentury, and then take the difference
# actually, the slopePerCentury will likely cancel out.....

# loop over longitude and latitude:
# find the linear best fit regression at every grid-point
for i in 1:dims[1]
    for j in 1:dims[2]
    #for j in 60:120
        #agrid[i,j],bgrid[i,j] = find_best_fit(sst_tr_mean[:],sst_summer[i,j,:])
        agrid[i,j],bgrid[i,j]   = find_best_fit(timeAxis[:],sst_summer[i,j,:])
        a_full[i,j],b_full[i,j] = find_best_fit(timeAxisFull[:],sst_var[i,j,:])
        #a_nino[i,j],b_nino[i,j] = find_best_fit(timeAxisFull[:],sst_var[i,j,:])
    end
end

# timeAxis should be a monthly time series, but with only 6 months per year.
# find slope of tropical mean sst, in celcius per year
a,b = find_best_fit(timeAxis,sst_tr_mean)
slopePerCentury = a*6*100
relSST = agrid*6*100 .- slopePerCentury   

# computing a 'relative' enso sst probably doesn't make sense.   what should the 
# independent variable be?
#relNinoSST = a_nino*600 .- slopePerCentury
#relNinaSST = a_nina*600 .- slopePerCentury

# compute the composite ENSO SST field:
#rel_ENSO_comp = relNinoSST .- relNinaSST
enso_high_mn  = mean(enso_high, dims=3);
enso_low_mn   = mean(enso_low, dims=3);
sst_ENSO_comp = enso_high_mn .- enso_low_mn

println("**************************************")
print("max value is: ",maximum(skipmissing(relSST)))
print("min value is: ",minimum(skipmissing(relSST)))
println("**************************************")
print("slope Per Century is: ",slopePerCentury)
println("**************************************")

#--------------------------------------------

## create arrays for the composites of the PI
#dims3     = size(vmax)
#pi_high  = Array{Union{Missing, Float64}, 3}(undef, dims3[1], dims3[2], numfields)
#pi_low   = Array{Union{Missing, Float64}, 3}(undef, dims3[1], dims3[2], numfields)
#
#for i in 1:numfields
#  pi_low[:,:,i]  = vmax[:,:,low[i]]    # La Nina months
#  pi_high[:,:,i] = vmax[:,:,high[i]]   # El Nino months
#end
#
#blackbird = pi_high .- pi_low;
#PI_comp_mn = mean(blackbird, dims=3)

#----------------------------------------------------------------------------------
## work on plots

levs = range(-10, 10, length = 21)
titsuf = "Nino - Nina"
#tit = tag * titsuf
tit = "MPI: " * titsuf
fig = Figure(;
    size = (800,900),
    )
    ax = GeoAxis(fig[2,1];
      xticks = -180:30:180,
      yticks = -90:30:90,
      limits=(-180,180,-40,40),
      #xlabel = "time (yr)",
      #ylabel = "RONI",
      title  = "Relative SST trend (K/Century)")
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
      title  = tag * " ENSO Threshold: " * thresh)
    lines!(C,ts_roni_sm, linestyle = :solid, linewidth=1.25)
    lines!(C,ensoPlotp, color = :black, linestyle = :solid, linewidth=1.25)
    lines!(C,ensoPlotn, color = :black,  linestyle = :solid, linewidth=1.25)
#plot!(C,ts_roni_sm) # plots the time series using dots for the data points. 
    plot!(C,forNew)
ax = GeoAxis(fig[3,1];
      xticks = -180:30:180,
      yticks = -90:30:90,
      limits=(-180,180,-40,40),
      title  = "Composite SST: Nino-Nina")
      bb = contourf!(ax, lon1d, lat1d, sst_ENSO_comp[:,:,1],
               levels = range(-2, 2, length = 21), 
               colormap = :vik,
               extendlow = :auto, extendhigh = :auto
          )
    lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=1.25)
    Colorbar(fig[3,2], bb)
ax = GeoAxis(fig[4,1];
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
    Colorbar(fig[4,2], bb)
fig


figname=tag*"_mpi_4pan.png"
#figname=tag*"_mpi_40x360_0p6K_150.png"
save(figname, fig)


##----------------------------------------------------------------------------
### save output to a netcdf file
#ds = NCDataset("test.nc","c")
#
#defDim(ds,"lon",144)
#defDim(ds,"lat",73)
#defDim(ds,"time",1)
#
#ds.attrib["title"] = "this is a test file"
#
#v1 = defVar(ds,"MPI",Float32,("lon","lat","time"))
#
#PI_comp_mn_nan = nomissing(PI_comp_mn,NaN)
#
#v1[:,:,:] = PI_comp_mn_nan
#
#v1.attrib["units"] = "meters per second"
#v1.attrib["commments"] = "this is the time avegage maximum potential intensity for TCs"
##v1.attrib["_FillValue"] = "Missing"
#
#close(ds)
