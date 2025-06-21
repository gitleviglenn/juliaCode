#-----------------------------------------------------------------------------
# mmm_analysis_cmip.jl
#
# read in SST (tos) and MPI(usually from scripts run on Casper)
#
# This script does a lot: 
#
# - compute a timeseries of RONI from the SST
# - compute the relative SST
# - select appropriate points for the positive and negative ENSO phases
# - plot the composite RH field for each model
# - plot the composite MPI field for each model
# - compute the trend line (slope) of RONI, for each model
# - compute the mean trend line for all models.   make plot.
#
# changes from tc_analysis_cmip.jl:  we now use the sst from within the 
# model_full_output_##.nc file, rather than a separate sst file.  as a result
# the input file to the function that calculates roni is not the full_output
# files, and that function now expects the sst variable name instead of the 
# tos variable name.  
# 
# to run from a terminal window:  
# >julia mmm_analysis_cmip.jl
#
# The work of Vecchi and Soden, 2007 focused on summer time, or NH TC season (June-November)
# months.   We can do that.   But if we want to look at ENSO based composites it makes more
# sense to start with the complete time series and then pick out the positive and negative
# phases.  
#
# ultimately it would be nice to output the MPI composite data as a netcdf file so that it can be 
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
using DataFrames

include("ensoFuncs.jl")
include("spec2relHumidity.jl") # output is relHcomp

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

# these values only work on a 360x180 grid
#tropS = 50
#tropN = 130
# these values only work on a 144x73 grid
tropS = 21
tropN = 53


timelen      = 1032
timelenH     = 516
timeAxis     = collect(1.:1:timelenH); 
timeAxisFull = collect(1.:1:timelen); 
# define an array that mirrors the time period of the data
C = collect(2015.083333:1/12:2101)
#-------------------------------------------------------------------------------------------------------------------
# specify file paths and names:

lpath="/Users/C823281551/data/cmip6/potInt/"


#file1   = lpath*"MPIESM/potInt/tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012.73x144.nc"
#file1b  = lpath*"MPIESM/potInt/tos_Omon_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501-210012_summer.73x144.nc"
#file2   = lpath*"MPIESM/potInt/MPI_ESM1_full_output.nc"

# MPI-ESM1-2-LR    First Model
#file1   = lpath*"MPI_ESM1_2_full_output_summer.nc"
#file2   = lpath*"MPI_ESM1_2_full_output_record.nc"

# MIROC6           Second Model
#file1   = lpath*"MIROC6_full_output_summer.nc"
#file2   = lpath*"MIROC6_full_output_record.nc"

# CESM2            Third Model                                      4.
#file1   = "/Users/C823281551/data/cmip6/CESM2/potInt/CESM_full_output.nc"
# ensemble member r10i1p1f1, computed on luft
#file2   = "/Users/C823281551/data/cmip6/CESM2/potInt/CESM_full_output.nc"
# ensemble member r4i1p1f1, computed on casper
#file2   = lpath*"CESM2_full_output.nc"

# Access-CM2       Fourth Model                                     1.
#file1   = lpath*"ACCESS_CM2_full_output_summer.nc"
#file2   = lpath*"ACCESS_CM2_full_output_record.nc"

# MRI-ESM2         Fifth Model
#file1   = lpath*"MPI_MRI-ESM2_full_output_summer.nc"
#file2   = lpath*"MPI_MRI-ESM2_full_output_record.nc"

# GFDL-ESM4        Sixth Model                                     10.
#file1   = lpath*"GFDL-ESM4_full_output_nohalo_summer.nc"
#file2   = lpath*"GFDL-ESM4_full_output_nohalo_record.nc"

# IPSL-CM6A        Seventch Model                                  9.
#file1   = lpath*"MPI_IPSL-CM6A_full_output_summer.nc"
#file2   = lpath*"MPI_IPSL-CM6A_full_output_record.nc"

# FGOALS-f3-L      Eighth Model                                    8.
#file2    = lpath*"FGOALS-f3-L_full_output.nc"

# MIROC6-ES2L      Ninth Model
#file2    = lpath*"MPI_MIROC-ES2L_full_output.nc"

# CanESM5          Tenth Model                                     3.
#file2    = lpath*"CanESM5_full_output_sethaloTOS.nc"

# CESM2-WACCM      Eleventh Model                                  6.
#file2    = lpath*"MPI_CESM2-WACCM_full_output.nc"

# CESM2-FV2        Twelfth Model                                   5.
#file2   = lpath*""

# EC-Earth3        Thirteenth Model                                7.
#file2   = lpath*"ECEarth3_full_output.nc"

# Access-ESMl-5    Fourteenth Model                                 2.
#file2   = lpath*"ACCESS_ESM15_full_output.nc"

# UKESM1-0-LL      Fifteenth Model
#file2   = lpath*"MPI_UKESM1-0-LL_full_output.nc"

# CNRM-ESM2-1      Sixteenth Model
#file2   =  lpath*"CNRM_ESM2_full_output.nc"

# HadGEM3-GC31-LL  Seventeenth Model

# GFDL-CM4                                                        11.                 

# NorESM2-LM --> still need to run analysis on casper

# FGOALS-g3

# AWI-CM-1-1-MR ??

# BCC-CSM2-MR ??

#-------------------------------------------------------------------------------------------------------------------
# tag to be used in plot title and figure name:

#------------------------------#
# First MODEL
#------------------------------#
file2   = lpath*"MPI_ESM1_2_full_output_record.nc"

tag = "MPI_ESM1_2"
#tag = "ACCESS_CM2"
#tag = "MIROC6"
#tag = "MRI_ESM2"
#tag = "GFDL_ESM4"
#tag = "MPI_IPSL"

## define parameters
## below are indices for a 360x40 degree resolution file from MPI-ESM: testing only
##lat1 = 1;lat2 = 40;lon34a = 10;lon34b = 61;lat34a = 15;lat34b = 25
#

# load the data
#data1 = NCDataset(file1)
#data1 = NCDataset(file1)
data2 = NCDataset(file2)
#data2b = NCDataset(file2b)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]
#lat1d = data1["lat"]
#lon1d = data1["lon"]

println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1")
#println("size of temp: ",size(temp))
#println("size of specific humidity: ",size(shum))
#shumtemp = shum[73,36,:,1]
#println("specific humidity: ",maximum(shumtemp))
println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1")

sst_var     = data2["sst"]
#sst_summer  = data1["sst"]
dims        = size(sst_var)
#dims2       = size(sst_summer)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);


#-------------------------------------------------------------------------------
# RONI
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
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
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

#RoniSlope1 = forNew

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

PI_comp_mn1 = PI_comp_mn
#-------------------------------------------------------------------------------
println(" working on relative SST calculation")

sst_nina = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
sst_nino = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
sst_trm_nina      = Array{Union{Missing, Float64}, 1}(undef, numfields)
sst_trm_nino      = Array{Union{Missing, Float64}, 1}(undef, numfields)

#lat1 = 50
#lat2 = 130

#latN = 120
#latS = 60

for i in 1:numfields
    #sst_tr_mean[i] = mean(skipmissing(sst_summer[:,tropS:tropN,i]))
    sst_trm_nina[i] = mean(skipmissing(sst_var[:,:,pad_ind_nina[i]]))
    sst_trm_nino[i] = mean(skipmissing(sst_var[:,:,pad_ind_nino[i]]))
end

for i in 1:numfields
  #enso_low[:,:,i]  = sst_var[:,:,pad_indices_Nina[i]]    # La Nina months
  #enso_high[:,:,i] = sst_var[:,:,pad_indices_Nino[i]]    # El Nino months
    sst_nino[:,:,i] = sst_var[:,:,pad_ind_nino[i]] .- sst_trm_nino[i]
    sst_nina[:,:,i] = sst_var[:,:,pad_ind_nina[i]] .- sst_trm_nina[i]
end   

#nino_composite = nino_mn - nina_mn
nino_composite  = sst_nino .- sst_nina
RelSST_comp_mn = mean(nino_composite, dims=3)

#RelSST_comp_mn1 = RelSST_comp_mn

## I think that I need to compute the relative sst change for Nino and Nina 
## separately, then subtract the slopePerCentury, and then take the difference
## actually, the slopePerCentury will likely cancel out.....
#
## loop over longitude and latitude:
## find the linear best fit regression at every grid-point
#for i in 1:dims[1]
#    for j in 1:dims[2]
#    #for j in 60:120
#        #agrid[i,j],bgrid[i,j] = find_best_fit(sst_tr_mean[:],sst_summer[i,j,:])
#        agrid[i,j],bgrid[i,j]   = find_best_fit(timeAxis[:],sst_summer[i,j,:])
#        a_full[i,j],b_full[i,j] = find_best_fit(timeAxisFull[:],sst_var[i,j,:])
#        #a_nino[i,j],b_nino[i,j] = find_best_fit(timeAxisFull[:],sst_var[i,j,:])
#    end
#end
#
## timeAxis should be a monthly time series, but with only 6 months per year.
## find slope of tropical mean sst, in celcius per year
#a,b = find_best_fit(timeAxis,sst_tr_mean)
#slopePerCentury = a*6*100
#relSST = agrid*6*100 .- slopePerCentury   

## compute the composite ENSO SST field:
#enso_high_mn  = mean(enso_high, dims=3);
#enso_low_mn   = mean(enso_low, dims=3);
#sst_ENSO_comp = enso_high_mn .- enso_low_mn

#println("**************************************")
#print("max value is: ",maximum(skipmissing(relSST)))
#print("min value is: ",minimum(skipmissing(relSST)))
#println("**************************************")
#print("slope Per Century is: ",slopePerCentury)
#println("**************************************")


#----------------------------------------------------------------------------------
# work on computing the relative humidity
println(" working on relative humidity calculation")

#nlats = 17
println("levels are: ",lev[:])
println("level 4 is: ",lev[4])

# 17 lat points, how many time points? 
#shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp_trH  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
temp_trL  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trH     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trL     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
num       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
den       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
svp       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhH       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhL       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

println("size of temp is: ",size(temp))
println("size of temp_tr is: ",size(temp_trH))
for i in 1:numfields
  #ji = 1
  #for j in 1:dims[2]
    temp_trH[:,:,i] = temp[:,:,4,pad_ind_nino[i]] .- 273.15 # convert to C   
    q_trH[:,:,i]    = shum[:,:,4,pad_ind_nino[i]]   
    temp_trL[:,:,i] = temp[:,:,4,pad_ind_nina[i]] .- 273.15 # convert to C 
    q_trL[:,:,i]    = shum[:,:,4,pad_ind_nina[i]]   
    #ji = ji + 1
  #end
end   

#println("temperature is (should be celcius): ",temp_trH[:,10,1])
println("specific humidity is (should be kg/kg): ",q_trH[:,10,1])

tempT = temp_trH
qT    = q_trH
num = exp.(34.494 .- 4924.99./(tempT .+ 237.1))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
#
rhH  = num./den

tempT = temp_trL
qT    = q_trL
num = exp.(34.494 .- (4924.99./(tempT .+ 237.1)))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
#
rhL  = num./den

# compute the composite relative humidity field
rh_comp = rhH .- rhL

# these magnitudes don't make sense to me at all...
RelRH_comp_mn = 1000000 .* mean(rh_comp, dims=4)

## test against values from Huang for saturation pressure for 
## temperatures below freezing: 
#tempTa = [-100, -80, -60, -40, -20, 0]
#num1 = exp.(43.494 .- (6545.8./(tempTa .+ 278)))
#den1 = (tempTa .+ 868).^2
#svp1 = num1./den1
#
#println("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaahhhhhhhhh")
#println(svp1[:])
#println("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaahhhhhhhhh")

#----------------------------------------------------------------------------------
## work on plots

levs = range(-2, 2, length = 21) 
tit="CMIP Relative ENSO composite SST (K)"
fig = fig_anom_plot(RelSST_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_SST.png"
save(figname, fig)

levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)

#forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

Slope1     = forest[1]
RoniSlope1 = forNew
RelSST_comp_mn1 = RelSST_comp_mn
RelRH_comp_mn1 = RelRH_comp_mn
PI_comp_mn_1 = PI_comp_mn

#println("relative humidity is: ",RelRH_comp_mn[:,10,1])

levs = range(-10, 10, length = 21) 
tit="CMIP RH composite (%)"
fig = fig_anom_plot(RelRH_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_RH.png"
save(figname, fig)

fig = Figure(;
    size = (900,1100), # width by height
    )
    ax = Axis(fig[1,1];
      xlabel = "time (yr)",
      ylabel = "RONI",
      title  = tag * " ENSO Threshold: " * thresh,)
    lines!(C,ts_roni_sm, linestyle = :solid, linewidth=1.25)
#    lines!(C,ensoPlotp, color = :black, linestyle = :solid, linewidth=1.25)
#    lines!(C,ensoPlotn, color = :black,  linestyle = :solid, linewidth=1.25)
##plot!(C,ts_roni_sm) # plots the time series using dots for the data points. 
    plot!(C,forNew)
figname=tag*"_slope.png"
save(figname, fig)

#------------------------------------------------------------------------------------
#------------------------------#
# Second MODEL
#------------------------------#

#file1   = lpath*"MIROC6_full_output_summer.nc"
file2   = lpath*"MIROC6_full_output_record.nc"

tag = "MIROC6"

# load the data
#data1 = NCDataset(file1)
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]
#lat1d = data1["lat"]
#lon1d = data1["lon"]

sst_var     = data2["sst"]
#sst_summer  = data1["sst"]
dims        = size(sst_var)
#dims2       = size(sst_summer)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
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
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

Slope2     = forest[1]
RoniSlope2 = forNew

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

PI_comp_mn2 = PI_comp_mn

levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)

PI_comp_mn_2 = PI_comp_mn

println("Done with second model: ",tag)
println("**************************************")


PI_comp_mn_tot = 0.5 .* (PI_comp_mn_1 .+ PI_comp_mn_2)

levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn_tot[:,:,1],lon,lat,tit,levs)
figname=tag*"_test_MPI_mmm.png"

save(figname, fig)
#-------------------------------------------------------------------------------
#------------------------------#
# Third MODEL
#------------------------------#

tag = "CESM2"

# ensemble member r4i1p1f1, computed on casper
file2   = lpath*"CESM2_full_output.nc"
#file2   = lpath*"CESM2_full_output_mh0.nc"
# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]


sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
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
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

Slope3     = forest[1]
RoniSlope3 = forNew

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

PI_comp_mn3 = PI_comp_mn

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")

#nlats = 17
println("levels are: ",lev[:])
println("level 4 is: ",lev[4])

# 17 lat points, how many time points? 
#shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp_trH  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
temp_trL  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trH     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trL     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
num       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
den       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
svp       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhH       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhL       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

println("size of temp is: ",size(temp))
println("size of temp_tr is: ",size(temp_trH))
for i in 1:numfields
  #ji = 1
  #for j in 1:dims[2]
    temp_trH[:,:,i] = temp[:,:,4,pad_ind_nino[i]] .- 273.15 # convert to C   
    q_trH[:,:,i]    = shum[:,:,4,pad_ind_nino[i]]   
    temp_trL[:,:,i] = temp[:,:,4,pad_ind_nina[i]] .- 273.15 # convert to C 
    q_trL[:,:,i]    = shum[:,:,4,pad_ind_nina[i]]   
    #ji = ji + 1
  #end
end   

#println("temperature is (should be celcius): ",temp_trH[:,10,1])
println("specific humidity is (should be kg/kg): ",q_trH[:,10,1])

tempT = temp_trH
qT    = q_trH
num = exp.(34.494 .- 4924.99./(tempT .+ 237.1))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from high ENSO values
rhH  = num./den

tempT = temp_trL
qT    = q_trL
num = exp.(34.494 .- (4924.99./(tempT .+ 237.1)))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from low ENSO values
rhL  = num./den

# compute the composite relative humidity field
rh_comp = rhH .- rhL

# these magnitudes don't make sense to me at all...
RelRH_comp_mn = 1000000 .* mean(rh_comp, dims=4)

RelRH_comp_mn3 = RelRH_comp_mn
#------------------------------------------------------------------------------------
# figures
levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)

#levs = range(-10, 10, length = 21) 
#tit="CMIP MPI composite (m/s)"
#fig = fig_anom_plot(PI_comp_mn_tot[:,:,1],lon,lat,tit,levs)
#figname=tag*"_test_MPI_mmm.png"

levs = range(-10, 10, length = 21) 
tit="CMIP RH composite (%)"
fig = fig_anom_plot(RelRH_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_RH.png"
save(figname, fig)

println("Done with the Third model: ",tag)
println("**************************************")

#------------------------------------------------------------------------------------
#------------------------------#
# Fourth MODEL
#------------------------------#

tag = "Access-CM2"
file2   = lpath*"ACCESS_CM2_full_output_record.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

Slope4     = forest[1]
RoniSlope4 = forNew

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

PI_comp_mn4 = PI_comp_mn
#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")

#nlats = 17
println("levels are: ",lev[:])
println("level 4 is: ",lev[4])

# 17 lat points, how many time points? 
#shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp_trH  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
temp_trL  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trH     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trL     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
num       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
den       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
svp       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhH       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhL       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

println("size of temp is: ",size(temp))
println("size of temp_tr is: ",size(temp_trH))
for i in 1:numfields
  #ji = 1
  #for j in 1:dims[2]
    temp_trH[:,:,i] = temp[:,:,4,pad_ind_nino[i]] .- 273.15 # convert to C   
    q_trH[:,:,i]    = shum[:,:,4,pad_ind_nino[i]]   
    temp_trL[:,:,i] = temp[:,:,4,pad_ind_nina[i]] .- 273.15 # convert to C 
    q_trL[:,:,i]    = shum[:,:,4,pad_ind_nina[i]]   
    #ji = ji + 1
  #end
end   

#println("temperature is (should be celcius): ",temp_trH[:,10,1])
println("specific humidity is (should be kg/kg): ",q_trH[:,10,1])

tempT = temp_trH
qT    = q_trH
num = exp.(34.494 .- 4924.99./(tempT .+ 237.1))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from high ENSO values
rhH  = num./den

tempT = temp_trL
qT    = q_trL
num = exp.(34.494 .- (4924.99./(tempT .+ 237.1)))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from low ENSO values
rhL  = num./den

# compute the composite relative humidity field
rh_comp = rhH .- rhL

# these magnitudes don't make sense to me at all...
RelRH_comp_mn = 1000000 .* mean(rh_comp, dims=4)

RelRH_comp_mn4 = RelRH_comp_mn

#------------------------------------------------------------------------------------
# figures
levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)

levs = range(-10, 10, length = 21) 
tit="CMIP RH composite (%)"
fig = fig_anom_plot(RelRH_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_RH.png"
save(figname, fig)

println("Done with the Fourth model: ",tag)
println("**************************************")
#------------------------------------------------------------------------------------
#------------------------------#
# Fifth MODEL
#------------------------------#
tag = "MRI-ESM2"
file2   = lpath*"MPI_MRI-ESM2_full_output_record.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")

#nlats = 17
println("levels are: ",lev[:])
println("level 4 is: ",lev[4])

# 17 lat points, how many time points? 
#shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp_trH  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
temp_trL  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trH     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trL     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
num       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
den       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
svp       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhH       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhL       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

println("size of temp is: ",size(temp))
println("size of temp_tr is: ",size(temp_trH))
for i in 1:numfields
  #ji = 1
  #for j in 1:dims[2]
    temp_trH[:,:,i] = temp[:,:,4,pad_ind_nino[i]] .- 273.15 # convert to C   
    q_trH[:,:,i]    = shum[:,:,4,pad_ind_nino[i]]   
    temp_trL[:,:,i] = temp[:,:,4,pad_ind_nina[i]] .- 273.15 # convert to C 
    q_trL[:,:,i]    = shum[:,:,4,pad_ind_nina[i]]   
    #ji = ji + 1
  #end
end   

#println("temperature is (should be celcius): ",temp_trH[:,10,1])
println("specific humidity is (should be kg/kg): ",q_trH[:,10,1])

tempT = temp_trH
qT    = q_trH
num = exp.(34.494 .- 4924.99./(tempT .+ 237.1))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from high ENSO values
rhH  = num./den

tempT = temp_trL
qT    = q_trL
num = exp.(34.494 .- (4924.99./(tempT .+ 237.1)))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from low ENSO values
rhL  = num./den

# compute the composite relative humidity field
rh_comp = rhH .- rhL

# these magnitudes don't make sense to me at all...
RelRH_comp_mn = 1000000 .* mean(rh_comp, dims=4)

# save fields from this model:
Slope5     = forest[1]
RoniSlope5     = forNew
PI_comp_mn5    = PI_comp_mn
RelRH_comp_mn5 = RelRH_comp_mn

#------------------------------------------------------------------------------------
# figures
levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)

levs = range(-10, 10, length = 21) 
tit="CMIP RH composite (%)"
fig = fig_anom_plot(RelRH_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_RH.png"
save(figname, fig)

println("Done with the Fifth model: ",tag)
println("**************************************")
#------------------------------#
# Sixth MODEL
#------------------------------#
tag = "GFDL-ESM4"
file2   = lpath*"GFDL-ESM4_full_output_nohalo_record.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")

#nlats = 17
println("levels are: ",lev[:])
println("level 4 is: ",lev[4])

# 17 lat points, how many time points? 
#shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp_trH  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
temp_trL  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trH     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trL     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
num       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
den       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
svp       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhH       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhL       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

println("size of temp is: ",size(temp))
println("size of temp_tr is: ",size(temp_trH))
for i in 1:numfields
  #ji = 1
  #for j in 1:dims[2]
    temp_trH[:,:,i] = temp[:,:,4,pad_ind_nino[i]] .- 273.15 # convert to C   
    q_trH[:,:,i]    = shum[:,:,4,pad_ind_nino[i]]   
    temp_trL[:,:,i] = temp[:,:,4,pad_ind_nina[i]] .- 273.15 # convert to C 
    q_trL[:,:,i]    = shum[:,:,4,pad_ind_nina[i]]   
    #ji = ji + 1
  #end
end   

#println("temperature is (should be celcius): ",temp_trH[:,10,1])
println("specific humidity is (should be kg/kg): ",q_trH[:,10,1])

tempT = temp_trH
qT    = q_trH
num = exp.(34.494 .- 4924.99./(tempT .+ 237.1))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from high ENSO values
rhH  = num./den

tempT = temp_trL
qT    = q_trL
num = exp.(34.494 .- (4924.99./(tempT .+ 237.1)))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from low ENSO values
rhL  = num./den

# compute the composite relative humidity field
rh_comp = rhH .- rhL

# these magnitudes don't make sense to me at all...
RelRH_comp_mn = 1000000 .* mean(rh_comp, dims=4)

# save fields from this model:
Slope6     = forest[1]
RoniSlope6     = forNew
PI_comp_mn6    = PI_comp_mn
RelRH_comp_mn6 = RelRH_comp_mn


#------------------------------------------------------------------------------------
# figures
levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)

levs = range(-10, 10, length = 21) 
tit="CMIP RH composite (%)"
fig = fig_anom_plot(RelRH_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_RH.png"
save(figname, fig)

println("Done with the Sixth model: ",tag)
println("**************************************")

#------------------------------#
# Seventh MODEL
#------------------------------#

tag = "IPSL-CM6A"
# IPSL-CM6A        Seventch Model
file2   = lpath*"MPI_IPSL-CM6A_full_output_record.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")

#nlats = 17
println("levels are: ",lev[:])
println("level 4 is: ",lev[4])

# 17 lat points, how many time points? 
#shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp_trH  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
temp_trL  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trH     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trL     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
num       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
den       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
svp       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhH       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhL       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

println("size of temp is: ",size(temp))
println("size of temp_tr is: ",size(temp_trH))
for i in 1:numfields
  #ji = 1
  #for j in 1:dims[2]
    temp_trH[:,:,i] = temp[:,:,4,pad_ind_nino[i]] .- 273.15 # convert to C   
    q_trH[:,:,i]    = shum[:,:,4,pad_ind_nino[i]]   
    temp_trL[:,:,i] = temp[:,:,4,pad_ind_nina[i]] .- 273.15 # convert to C 
    q_trL[:,:,i]    = shum[:,:,4,pad_ind_nina[i]]   
    #ji = ji + 1
  #end
end   

#println("temperature is (should be celcius): ",temp_trH[:,10,1])
println("specific humidity is (should be kg/kg): ",q_trH[:,10,1])

tempT = temp_trH
qT    = q_trH
num = exp.(34.494 .- 4924.99./(tempT .+ 237.1))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from high ENSO values
rhH  = num./den

tempT = temp_trL
qT    = q_trL
num = exp.(34.494 .- (4924.99./(tempT .+ 237.1)))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from low ENSO values
rhL  = num./den

# compute the composite relative humidity field
rh_comp = rhH .- rhL

# these magnitudes don't make sense to me at all...
RelRH_comp_mn = 1000000 .* mean(rh_comp, dims=4)

# save fields from this model:
Slope7         = forest[1]
RoniSlope7     = forNew
PI_comp_mn7    = PI_comp_mn
RelRH_comp_mn7 = RelRH_comp_mn

#-------------------------------------------------------------------------------------------------------------------
# MIROC6-ES2L      Ninth Model
#file2    = lpath*"MPI_MIROC-ES2L_full_output.nc"

# CanESM5          Tenth Model
#file2    = lpath*"CanESM5_full_output_sethaloTOS.nc"

# CESM2-WACCM      Eleventh Model
#file2    = lpath*"MPI_CESM2-WACCM_full_output.nc"
#-------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# figures
levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)

levs = range(-10, 10, length = 21) 
tit="CMIP RH composite (%)"
fig = fig_anom_plot(RelRH_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_RH.png"
save(figname, fig)

println("Done with the Seventh model: ",tag)
println("**************************************")

#------------------------------------------------------------------------------------
#------------------------------#
# Eighth MODEL
#------------------------------#

tag = "FGOALS-f3-L"
file2    = lpath*"FGOALS-f3-L_full_output.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")

#nlats = 17
println("levels are: ",lev[:])
println("level 4 is: ",lev[4])

# 17 lat points, how many time points? 
#shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp_trH  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
temp_trL  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trH     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trL     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
num       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
den       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
svp       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhH       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhL       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

println("size of temp is: ",size(temp))
println("size of temp_tr is: ",size(temp_trH))
for i in 1:numfields
  #ji = 1
  #for j in 1:dims[2]
    temp_trH[:,:,i] = temp[:,:,4,pad_ind_nino[i]] .- 273.15 # convert to C   
    q_trH[:,:,i]    = shum[:,:,4,pad_ind_nino[i]]   
    temp_trL[:,:,i] = temp[:,:,4,pad_ind_nina[i]] .- 273.15 # convert to C 
    q_trL[:,:,i]    = shum[:,:,4,pad_ind_nina[i]]   
    #ji = ji + 1
  #end
end   

#println("temperature is (should be celcius): ",temp_trH[:,10,1])
println("specific humidity is (should be kg/kg): ",q_trH[:,10,1])

tempT = temp_trH
qT    = q_trH
num = exp.(34.494 .- 4924.99./(tempT .+ 237.1))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from high ENSO values
rhH  = num./den

tempT = temp_trL
qT    = q_trL
num = exp.(34.494 .- (4924.99./(tempT .+ 237.1)))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from low ENSO values
rhL  = num./den

# compute the composite relative humidity field
rh_comp = rhH .- rhL

# these magnitudes don't make sense to me at all...
RelRH_comp_mn = 1000000 .* mean(rh_comp, dims=4)

# save fields from this model:
Slope8         = forest[1]
roni_ts8       = ts_roni_sm
RoniSlope8     = forNew
PI_comp_mn8    = PI_comp_mn
RelRH_comp_mn8 = RelRH_comp_mn

levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)


println("Done with the Eighth model: ",tag)
println("**************************************")

#------------------------------------------------------------------------------------
#------------------------------#
# Ninth MODEL
#------------------------------#
tag = "MIROC6-ES2L"
file2    = lpath*"MPI_MIROC-ES2L_full_output.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")

#nlats = 17
println("levels are: ",lev[:])
println("level 4 is: ",lev[4])

# 17 lat points, how many time points? 
#shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp_trH  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
temp_trL  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trH     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trL     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
num       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
den       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
svp       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhH       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhL       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

println("size of temp is: ",size(temp))
println("size of temp_tr is: ",size(temp_trH))
for i in 1:numfields
  #ji = 1
  #for j in 1:dims[2]
    temp_trH[:,:,i] = temp[:,:,4,pad_ind_nino[i]] .- 273.15 # convert to C   
    q_trH[:,:,i]    = shum[:,:,4,pad_ind_nino[i]]   
    temp_trL[:,:,i] = temp[:,:,4,pad_ind_nina[i]] .- 273.15 # convert to C 
    q_trL[:,:,i]    = shum[:,:,4,pad_ind_nina[i]]   
    #ji = ji + 1
  #end
end   

#println("temperature is (should be celcius): ",temp_trH[:,10,1])
println("specific humidity is (should be kg/kg): ",q_trH[:,10,1])

tempT = temp_trH
qT    = q_trH
num = exp.(34.494 .- 4924.99./(tempT .+ 237.1))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from high ENSO values
rhH  = num./den

tempT = temp_trL
qT    = q_trL
num = exp.(34.494 .- (4924.99./(tempT .+ 237.1)))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from low ENSO values
rhL  = num./den

# compute the composite relative humidity field
rh_comp = rhH .- rhL

# these magnitudes don't make sense to me at all...
RelRH_comp_mn = 1000000 .* mean(rh_comp, dims=4)

# save fields from this model:
Slope9         = forest[1]
roni_ts9       = ts_roni_sm
RoniSlope9     = forNew
PI_comp_mn9    = PI_comp_mn
RelRH_comp_mn9 = RelRH_comp_mn

levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)


println("Done with the Ninth model: ",tag)
println("**************************************")

#-------------------------------------------------------------------------------------------------------------------

#------------------------------#
# Tenth MODEL
#------------------------------#
tag = "CanESM5"       #   Tenth Model
file2    = lpath*"CanESM5_full_output_sethaloTOS.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")

#nlats = 17
println("levels are: ",lev[:])
println("level 4 is: ",lev[4])

# 17 lat points, how many time points? 
#shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp_trH  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
temp_trL  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trH     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trL     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
num       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
den       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
svp       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhH       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhL       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

println("size of temp is: ",size(temp))
println("size of temp_tr is: ",size(temp_trH))
for i in 1:numfields
  #ji = 1
  #for j in 1:dims[2]
    temp_trH[:,:,i] = temp[:,:,4,pad_ind_nino[i]] .- 273.15 # convert to C   
    q_trH[:,:,i]    = shum[:,:,4,pad_ind_nino[i]]   
    temp_trL[:,:,i] = temp[:,:,4,pad_ind_nina[i]] .- 273.15 # convert to C 
    q_trL[:,:,i]    = shum[:,:,4,pad_ind_nina[i]]   
    #ji = ji + 1
  #end
end   

#println("temperature is (should be celcius): ",temp_trH[:,10,1])
println("specific humidity is (should be kg/kg): ",q_trH[:,10,1])

tempT = temp_trH
qT    = q_trH
num = exp.(34.494 .- 4924.99./(tempT .+ 237.1))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from high ENSO values
rhH  = num./den

tempT = temp_trL
qT    = q_trL
num = exp.(34.494 .- (4924.99./(tempT .+ 237.1)))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from low ENSO values
rhL  = num./den

# compute the composite relative humidity field
rh_comp = rhH .- rhL

# these magnitudes don't make sense to me at all...
RelRH_comp_mn = 1000000 .* mean(rh_comp, dims=4)

# save fields from this model:
Slope10         = forest[1]
roni_ts10       = ts_roni_sm
RoniSlope10     = forNew
PI_comp_mn10    = PI_comp_mn
RelRH_comp_mn10 = RelRH_comp_mn

levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)

println("Done with the Tenth model: ",tag)
println("**************************************")

#------------------------------#
# Eleventh MODEL
#------------------------------#
tag = "CESM2-WACCM"  #    Eleventh Model
file2    = lpath*"MPI_CESM2-WACCM_full_output.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")

#nlats = 17
println("levels are: ",lev[:])
println("level 4 is: ",lev[4])

# 17 lat points, how many time points? 
#shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp_trH  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
temp_trL  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trH     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
q_trL     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
num       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
den       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
svp       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhH       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
rhL       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

println("size of temp is: ",size(temp))
println("size of temp_tr is: ",size(temp_trH))
for i in 1:numfields
  #ji = 1
  #for j in 1:dims[2]
    temp_trH[:,:,i] = temp[:,:,4,pad_ind_nino[i]] .- 273.15 # convert to C   
    q_trH[:,:,i]    = shum[:,:,4,pad_ind_nino[i]]   
    temp_trL[:,:,i] = temp[:,:,4,pad_ind_nina[i]] .- 273.15 # convert to C 
    q_trL[:,:,i]    = shum[:,:,4,pad_ind_nina[i]]   
    #ji = ji + 1
  #end
end   

#println("temperature is (should be celcius): ",temp_trH[:,10,1])
println("specific humidity is (should be kg/kg): ",q_trH[:,10,1])

tempT = temp_trH
qT    = q_trH
num = exp.(34.494 .- 4924.99./(tempT .+ 237.1))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from high ENSO values
rhH  = num./den

tempT = temp_trL
qT    = q_trL
num = exp.(34.494 .- (4924.99./(tempT .+ 237.1)))
den = (tempT .+ 105).^1.57
svp = num./den
# reuse the num and den arrays:
num = qT
den = (1. .- qT).*svp
# relative humidity from low ENSO values
rhL  = num./den

# compute the composite relative humidity field
rh_comp = rhH .- rhL

# these magnitudes don't make sense to me at all...
RelRH_comp_mn = 1000000 .* mean(rh_comp, dims=4)

# save fields from this model:
Slope11         = forest[1]
roni_ts11       = ts_roni_sm
RoniSlope11     = forNew
PI_comp_mn11    = PI_comp_mn
RelRH_comp_mn11 = RelRH_comp_mn

levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)

println("Done with the Eleventh model: ",tag)
println("**************************************")

#------------------------------#
# Twelfth MODEL
#------------------------------#

# CESM2-FV2
tag = "CESM2-FV2"  #    Twelfth Model
file2    = lpath*"MPI_CESM2-FV2_full_output.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")
spec2relComp(lpath,file2,pad_ind_nino,pad_ind_nina)
RelRH_comp_mn = relHcomp

#------------------------------------------------------------------------------------
# save fields from this model:
Slope12         = forest[1]
roni_ts12       = ts_roni_sm
RoniSlope12     = forNew
PI_comp_mn12    = PI_comp_mn
RelRH_comp_mn12 = RelRH_comp_mn

println("Done with the Twelfth model: ",tag)
println("**************************************")

#------------------------------#
# Thirteenth MODEL
#------------------------------#

tag = "EC-Earth3"  #    Thirteenth Model
file2   = lpath*"ECEarth3_full_output.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")
spec2relComp(lpath,file2,pad_ind_nino,pad_ind_nina)
RelRH_comp_mn = relHcomp

#------------------------------------------------------------------------------------
# save fields from this model:
Slope13         = forest[1]
roni_ts13       = ts_roni_sm
RoniSlope13     = forNew
PI_comp_mn13    = PI_comp_mn
RelRH_comp_mn13 = RelRH_comp_mn

println("Done with the Thirteenth model: ",tag)
println("**************************************")

#------------------------------#
# Fourteenth MODEL
#------------------------------#

tag     = "Access-ESMl-5"     #  Fourteenth Model
file2   = lpath*"ACCESS_ESM15_full_output.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")
spec2relComp(lpath,file2,pad_ind_nino,pad_ind_nina)
RelRH_comp_mn = relHcomp

#------------------------------------------------------------------------------------
# save fields from this model:
Slope14         = forest[1]
roni_ts14       = ts_roni_sm
RoniSlope14     = forNew
PI_comp_mn14    = PI_comp_mn
RelRH_comp_mn14 = RelRH_comp_mn

println("Done with the Fourteenth model: ",tag)
println("**************************************")

#------------------------------#
# Fifteenth MODEL
#------------------------------#

tag = "UKESM1-0-LL"   #    Fifteenth Model
file2   = lpath*"MPI_UKESM1-0-LL_full_output.nc"

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")
spec2relComp(lpath,file2,pad_ind_nino,pad_ind_nina)
RelRH_comp_mn = relHcomp

#------------------------------------------------------------------------------------
# save fields from this model:
Slope15         = forest[1]
roni_ts15       = ts_roni_sm
RoniSlope15     = forNew
PI_comp_mn15    = PI_comp_mn
RelRH_comp_mn15 = RelRH_comp_mn

println("Done with the Fifteenth model: ",tag)
println("**************************************")



#------------------------------#
# Sixteenth MODEL
#------------------------------#

#------------------------------------------------------------------------------------
tag = "CNRM-ESM2-1"
file2   =  lpath*"CNRM_ESM2_full_output.nc"
#------------------------------------------------------------------------------------

# load the data
data2 = NCDataset(file2)
vmax  = data2["vmax"]
shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
temp  = data2["t"] # air temperature (K)
lat   = data2["lat"]
lon   = data2["lon"]
lev   = data2["plev"]

sst_var     = data2["sst"]
dims        = size(sst_var)

#sst_tr_mean = Array{Union{Missing, Float64}, 1}(undef, dims2[3])
agrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
a_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
b_full      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pi_high     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
pi_low      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_high   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);
enso_low    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields);

#-------------------------------------------------------------------------------
# RONI:
inpFile = file2 
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(" working on MPI calculation")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")

# call function to process cmip data and compute RONI:
calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn

# choose to use, or not, a running mean filter to smooth the timeseries.
#smooth_12_ts(ts_roni,timelen)
#ts_roni_sm = ts_12_sm
ts_roni_sm = ts_roni

# calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

println("**************************************")
println("slope of time series is: ",forest[1])
println("**************************************")

## This method identifies a number of ENSO events and then the function 
## padding selects 6 indices for each event.  
local_max_indices = find_local_maxima(ts_roni_sm,ensoDef,numOfevents)
local_min_indices = find_local_minima(ts_roni_sm,-ensoDef,numOfevents)

pad_ind_nino = padding(local_max_indices)
pad_ind_nina = padding(local_min_indices)
println("pad_indices_Nina: ",pad_ind_nina[:])

for i in 1:numfields
  pi_low[:,:,i]  = vmax[:,:,pad_ind_nina[i]]    # La Nina months
  pi_high[:,:,i] = vmax[:,:,pad_ind_nino[i]]   # El Nino months
end   

# compute the composite MPI field: 
pi_high_mn = mean(pi_high, dims=3);
pi_low_mn  = mean(pi_low, dims=3);
PI_comp_mn = pi_high_mn .- pi_low_mn;

#------------------------------------------------------------------------------------
# relative humidity
println(" working on relative humidity calculation")
spec2relComp(lpath,file2,pad_ind_nino,pad_ind_nina)
RelRH_comp_mn = relHcomp

#------------------------------------------------------------------------------------
# save fields from this model:
Slope16         = forest[1]
roni_ts16       = ts_roni_sm
RoniSlope16     = forNew
PI_comp_mn16    = PI_comp_mn
RelRH_comp_mn16 = RelRH_comp_mn


#------------------------------------------------------------------------------------

println("Done with the Sixteenth model: ",tag)
println("**************************************")

## test roni calcluation with different input tos file: 
#file2 = "/Users/C823281551/data/cmip6/potInt/CNRM_ESM2_tos_73x144_pp.nc"
#
#inpFile = file2 
#println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
#println(" working on MPI calculation")
#println(inpFile)
#println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
#
## call function to process cmip data and compute RONI:
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
##prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#ts_roni = ts_rmn
#
## choose to use, or not, a running mean filter to smooth the timeseries.
##smooth_12_ts(ts_roni,timelen)
##ts_roni_sm = ts_12_sm
#ts_roni_sm = ts_roni
#
## calculate the best fit line
#bb = collect(0:1031);
#ensoPlotp = zeros(1032) .+ ensoDef
#ensoPlotn = zeros(1032) .- ensoDef
#forest = find_best_fit(bb,ts_roni_sm)
#forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
## to save for each model
#
#Slope17         = forest[1]
#roni_ts17       = ts_roni_sm

#------------------------------#
# Seventeenth MODEL
#------------------------------#

#------------------------------------------------------------------------------------
tag = "HadGEM3-GC31-LL"
file2   = "/Users/C823281551/data/cmip6/HadGEM3/tos_73x144.nc" 
inpFile = file2
## call function to process cmip data and compute RONI:
#calc_roni_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
##prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ts_roni = ts_rmn
ts_roni_sm = ts_roni

## calculate the best fit line
bb = collect(0:1031);
ensoPlotp = zeros(1032) .+ ensoDef
ensoPlotn = zeros(1032) .- ensoDef
forest = find_best_fit(bb,ts_roni_sm)
forNew = collect(bb).*forest[1] .+ forest[2]; # this is the slope that we want 
# to save for each model

Slope17         = forest[1]
roni_ts17       = ts_roni_sm
#------------------------------------------------------------------------------------

println("Done with the Seventeenth model: ",tag)
println("**************************************")

#------------------------------#
# Eighthteenth MODEL
#------------------------------#



#------------------------------#
#------------------------------#
#------------------------------#

println("**************************************")
println("**************************************")
println("slopes from each of the models: ")
println("slope = ",Slope1)
println("slope = ",Slope2)
println("slope = ",Slope3)
println("slope = ",Slope4)
println("slope = ",Slope5)
println("slope = ",Slope6)
println("slope = ",Slope7)
println("slope = ",Slope8)
println("slope = ",Slope9)
println("slope = ",Slope10)
println("slope = ",Slope11)
println("slope = ",Slope12)
println("slope = ",Slope13)
println("slope = ",Slope14)
println("slope = ",Slope15)
println("slope = ",Slope16)
println("slope = ",Slope17)
println("**************************************")

total = Slope1 + Slope2 + Slope3 +  Slope4 +  Slope5 +  Slope6 +  Slope7 +  Slope8 +  Slope9 +  Slope10 +  Slope11 +  Slope12 +  Slope13 +  Slope14 +  Slope15 + Slope16 + Slope17

mn_slope = total/17
println("mean slope = ",mn_slope)
println("**************************************")
println("**************************************")

#-------------------------------------------------------------------------------------------------------------------

levs = range(-10, 10, length = 21) 
tit="CMIP MPI composite (m/s)"
fig = fig_anom_plot(PI_comp_mn[:,:,1],lon,lat,tit,levs)
figname=tag*"_MPI.png"
save(figname, fig)

# figure showing the time series and trend of RONI:
fig = Figure(;
    size = (900,1100), # width by height
    )
    ax = Axis(fig[1,1];
      xlabel = "time (yr)",
      ylabel = "RONI",
      title  = "RONI ENSO timeseries slopes: " * thresh,)
    lines!(C,roni_ts11, linestyle = :solid, linewidth=1.25)
    #lines!(C,ts_roni_sm, linestyle = :solid, linewidth=1.25)
#    lines!(C,ensoPlotp, color = :black, linestyle = :solid, linewidth=1.25)
#    lines!(C,ensoPlotn, color = :black,  linestyle = :solid, linewidth=1.25)
##plot!(C,ts_roni_sm) # plots the time series using dots for the data points. 
    lines!(C,RoniSlope1, linestyle = :solid, linewidth=1.5)
    lines!(C,RoniSlope2, linewidth=1.5)
    lines!(C,RoniSlope3, linewidth=1.5)
    lines!(C,RoniSlope4, linewidth=1.5)
    lines!(C,RoniSlope5, linewidth=1.5)
    lines!(C,RoniSlope6, linewidth=1.5)
    lines!(C,RoniSlope7, linewidth=1.5)
    lines!(C,roni_ts8, linestyle = :solid, linewidth=1.25)
    lines!(C,RoniSlope8, linewidth=1.5)
    lines!(C,roni_ts9, linestyle = :solid, linewidth=1.25)
    lines!(C,RoniSlope9, linewidth=1.5)
    lines!(C,roni_ts10, linestyle = :solid, linewidth=1.25)
    lines!(C,RoniSlope10, linewidth=1.5)
    lines!(C,RoniSlope11, linewidth=1.5)
figname="RONI_ts_slopes.png"
save(figname, fig)

# figure showing the time series and trend of RONI:
fig = Figure(;
    size = (900,1100), # width by height
    )
    ax = Axis(fig[1,1];
      xlabel = "time (yr)",
      ylabel = "RONI",
      title  = "RONI ENSO timeseries slopes: " * thresh,)
    lines!(C,roni_ts11, linestyle = :solid, linewidth=1.25)
    #lines!(C,ts_roni_sm, linestyle = :solid, linewidth=1.25)
#    lines!(C,ensoPlotp, color = :black, linestyle = :solid, linewidth=1.25)
#    lines!(C,ensoPlotn, color = :black,  linestyle = :solid, linewidth=1.25)
##plot!(C,ts_roni_sm) # plots the time series using dots for the data points. 
    lines!(C,roni_ts16, linestyle = :solid, linewidth=1.25)
    lines!(C,RoniSlope9, linewidth=1.5)
    lines!(C,roni_ts17, linestyle = :solid, linewidth=1.25)
figname="RONI_test_ts_slopes.png"
save(figname, fig)


