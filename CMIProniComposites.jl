#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CMIProniComposites.jl
#
# ultimate goal:  produce a field of RH that show the difference of El Nino and La Nina 
# states for a particular CMIP6 model
#
# create a timeseries of RONI, set a threshold that defines the positive and negative 
# phase of ENSO, and then based on that timeseries, grab the appropriate times in the 
# RH fields that correspond to either positive or negative phases of ENSO
#
# ensoFuncs.jl contains many of the functions used in this script and needs to be run 
# before this script will work.  
#
# both SST (tos) and RH (hur) fields need to be read into script
#
# to do:  1.  get rid of the asinine white column at 0 degrees longitude.  
#         2.  figure out how to shift the center longitude 
#
# Major Problem with MPI tos:  the prepare_cmip_ts function looks like the indices used 
#                              to define the nino3.4 region are model specific.  now that
#                              I am using a file with the entire globe, the box used is 
#                              almost certainly wrong
#
# this should be useful for creating plots with EArth data: 
# https://forem.julialang.org/gage/creating-nice-simple-geo-plots-from-scratch-using-makiejl-301
# 
# levi silvers                                                                nov 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using DataFrames
using CairoMakie
using GeoMakie
using Statistics 
using NCDatasets

include("ensoFuncs.jl")

# calculate RONI
# based on a threshold of RONI determine ENSO events
# grab RH fields when ENSO events are occuring
# plot the composite RH fields

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# beginning of work with GCM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path="/Users/C823281551/data/"

# what is the best way to process each model?   something like cases that change the paths?  
# simply changing the model name will not be good enough


modelp="CESM2"
filehur  = path*"cmip6/CESM2/hur_Amon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215_regrid.nc"
filetos  = path*"cmip6/CESM2/tos_Omon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215_regrid2.nc" 
tag = "CESM2"

#modelp="MPI-ESM1-2-LR"
###filehur  = path*"cmip6/MPIESM/hur_Amon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_regridded.nc"
#filehur  = path*"cmip6/MPIESM/hur_Amon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_regridded4.nc"
#filetos  = path*"cmip6/MPIESM/tos_Omon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_pm40b.nc" 
##filetos  = path*"cmip6/MPIESM/tos_Omon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_latlon.nc" 
##filetos = path*"MPI-ESM1/tos_Omon_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_20150116-21001216_latlon.nc"
#tag = "MPI"

#filehur  = path*"cmip6/CNRMESM2/hur_Amon_CNRM-ESM2-1_ssp585_r1i1p1f2_gr_20150116-21001216_regrid.nc"
#filetos  = path*"cmip6/CNRMESM2/tos_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_20150116-21001216_regrid2.nc" 
#tag = "CNRMESM2"

#filehur  = path*"cmip6/HadGEM3/hur_Amon_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_20150116-21001216_regridFull.nc"
#filetos  = path*"cmip6/HadGEM3/tos_Omon_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_20150116-21001216_regridFull.nc" 
#tag = "HadGEM3"

#filehur  = path*"cmip6/ACCESS/hur_Amon_ACCESS-CM2_ssp585_r1i1p1f1_gn_20150116-21001216_regrid.nc"
#filetos  = path*"cmip6/ACCESS/tos_Omon_ACCESS-CM2_ssp585_r1i1p1f1_gn_20150116-21001216_regrid.nc" 
#tag = "ACCESS"

###FGOALS
#filehur = path*"cmip6/FGOALS/hur_Amon_FGOALS-g3_ssp585_r1i1p1f1_gn_20150116-21001216_regrid.nc"
#filetos = path*"cmip6/FGOALS/tos_Omon_FGOALS-g3_ssp585_r1i1p1f1_gn_20150116-21001216_regrid.nc"
#tag = "FGOALS"

#file1b  = path*"cmip6/CESM2/tos_Omon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215.nc" 

#filehur  = path*"cmip6/MPIESM/hur_Amon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_regridded.nc"
#file1b  = path*"cmip6/MPIESM/tos_Omon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_latlon.nc" 

#filetos  = path*"cmip6/MPIESM/tos_Omon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_latlon.nc" 

data   = NCDataset(filehur)
data1  = NCDataset(filetos)

data.attrib

rh   = data["hur"] # hur(time, plev, lat, lon)
sst  = data1["tos"] # hur(time, plev, lat, lon)
latb = data1["lat"]
lat  = data["lat"]
lon  = data["lon"]
tme  = data["time"]
lev  = data["plev"]

#tim = 5 # we need to grab rh values at the times selected by check_thresh
#rh1 = rh[:,:,1,tim]
#rh2 = rh[:,:,2,tim]

#scenario timeseries
timelen2=1032
inpFile = filetos
#println("~~~~~~~~~~~~~~~~~~file 6b~~~~~~~~~~~~~~~~~~~~~~")
#println(inpFile)  # HadGEM3
#println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
lat1   = 71
lat2   = 110
lat34a = 85 
lat34b = 96 
lon34a = 10
lon34b = 61
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
# MPI
# these values work for a file that extends between +/-20 degrees
#lat1 = 1
#lat2 = 40
#lon34a = 10
#lon34b = 61
#lat34a = 15
#lat34b = 25

## works for MPI for +/-40 degrees 
#lat1 = 21 
#lat2 = 60
#lon34a = 10
#lon34b = 61
#lat34a = 35
#lat34b = 46

## FGOALS
#lat1 = 71 
#lat2 = 110
#lon34a = 10
#lon34b = 61
#lat34a = 85
#lat34b = 96

prepare_cmip_ts(inpFile,timelen2,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba1b = ts_rmn

# smooth time series to be checked for threshold
smooth_12_ts(ba1b,timelen2)
ba1b_sm1 = ts_12_sm

# HadGEM3 has an enso time series that has a strong negative trend...
# ensoDef = 1.0; HadGEM3
ensoDef = 1.7
# MPI ensoDef = 1.7
thshd = ensoDef

# i think the only reason we need to pass in a file is to get the 
# time dimension from it.   to speed things up we could pass that
# in instead of dealing with the file within check_thresh()
check_thresh_high(inpFile, ba1b_sm1, thshd) # output is 'high'
thshd = -ensoDef
#thshd = -2.0 # FGOALS
check_thresh_low(inpFile, ba1b_sm1, thshd)  # output is 'low'

# read in RH file
ds = NCDataset(filehur)
ds.attrib
nctime = ds["time"]
println("~~~~~~~~~Silence~~~~~~~~~~~~~~~~~~")
println(nctime[high[1:8]])
println("~~~~~~~~~Golden~~~~~~~~~~~~~~~~~~~")

#rhhigh = rh[:,:,1,high]

dims      = size(rh)
dims2     = size(sst)
#numfields = 8 # for FGOALS this was 50...  45 seems to work well for MPIESM... 
numfields = 40
# thshd = 1.5 and numfields = 8 works best for FGOALS...
rh_high  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
rh_low   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
sst_high = Array{Union{Missing, Float64}, 3}(undef, dims2[1], dims2[2], numfields)
sst_low  = Array{Union{Missing, Float64}, 3}(undef, dims2[1], dims2[2], numfields)

#-------------------------------------------
# quick test: 
#
# restults can be plotted in the REPL as: 
# test_rh = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 2)
# test_rh[:,:,1] = rh[:,:,1,high[1]]
# test_rh[:,:,2] = rh[:,:,1,high[2]]
# contourf(test_rh[:,:,1])
#-------------------------------------------

# RH is computed as the geometric mean between 2 vertical levels
endi = numfields
for i in 1:endi
  rh_low[:,:,i]   = 0.5(rh[:,:,1,low[i]].+rh[:,:,2,low[i]])
  sst_low[:,:,i]  = sst[:,:,low[i]]
  #rh_high[:,:,i] = rh[:,:,1,high[i]]
  sst_high[:,:,i]  = sst[:,:,high[i]]
  rh_high[:,:,i]   = 0.5(rh[:,:,1,high[i]].+rh[:,:,2,high[i]])
end

# average in time
rh_high_mn  = mean(rh_high, dims = 3)
rh_low_mn   = mean(rh_low, dims = 3)
sst_high_mn = mean(sst_high, dims = 3)
sst_low_mn  = mean(sst_low, dims = 3)

# compute the difference between the high and low ENSO indices
rh_diff  = rh_high_mn - rh_low_mn
sst_diff = sst_high_mn - sst_low_mn

# get rh fields for positive enso phase (above thshd) 
# get rh fields for negative enso phase
# average the rh fields for each phase
# compute and plot the difference of the fields, El Nino - La Nina


function fig_1_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(600,300),
        )
    ax = Axis(f2[1,1];
        #xticks = -180:30:180, 
        xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             levels = range(0, 100, length = 10), 
             colormap = :Blues_8,
             #colormap = :navia,
             #colormap = :batlow,
             #colormap = :vik,
        )
        Colorbar(f2[1,2], bb)
    return f2
end
# define the limits of the plotted area
lon1=-180
lon2=180
lat1=-40
lat2=40
function fig_anom_plot(inpv,d1,d2,tit,levs)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        #backgroundcolor=:snow2,
        backgroundcolor=:white,
        size=(600,300),
        )
    ax = GeoAxis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        xlabel="longitude",
        ylabel="latitude",
        limits=(lon1,lon2,lat1,lat2),
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             levels = levs,
             #levels = range(-20, 20, length = 100),
             #colormap = :batlow,
             colormap = :vik, # was default for redish bluish
             #colormap = :BrBg, # better for RH  browns and greens
             #colormap = :roma,
             extendlow = :auto, extendhigh = :auto
        )
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[1,2], bb)
    return f2
end
function fig_2_plot(inpv,d1,d2,tit)
    f = Figure(;
        #figure_padding=(5,5,10,10),
        #backgroundcolor=:snow2,
        size=(600,300),
        )
    ax = Axis(f[1,1];
        #xticks = -180:30:180, 
        ##xticks = 0:30:360, 
        #yticks = -90:30:90,
        #xlabel="longitude",
        #ylabel="latitude",
        #limits=(lon1,lon2,lat1,lat2),
        #title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             levels = range(-20, 20, length = 100),
             #colormap = :batlow,
             colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        #lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f[1,2], bb)
    return f
end

#rh_2_plot = rh_low_mn[:,:,1]
#fig = fig_1_plot(rh_2_plot,lon,lat,"RH low")

#rh_2_plot = rh_high_mn[:,:,1]
#fig = fig_1_plot(rh_2_plot,lon,lat,"RH high")

# figure 1
rh_anom = rh_diff[:,:,1]
tit1="RH anomaly, "*tag*" ssp585"
fig1name=tag*"rh.png"
level1 = range(-20, 20, length = 21)
fig = fig_anom_plot(rh_anom,lon,lat,tit1,level1)
save(fig1name, fig)

# figure 2
sst_anom = sst_diff[:,:,1]
#contourf(sst_anom[:,:,1])
tit2="SST anomaly, "*tag*" ssp585"
level1 = range(-5, 5, length = 50)
#fig = fig_anom_plot(sst_anom,lon,latb,tit2,level1)
fig = fig_anom_plot(sst_anom,lon,latb,tit2,level1)
fig2name=tag*"sst.png"
save(fig2name, fig)



