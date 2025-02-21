#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotENSOcmipIndices.jl
#
# ensoFuncs.jl contains many of the functions used in this script and needs to be run 
# before this script will work.  
#
# - plot ENSO indices for individual cmip models
# - nino 3.4 region is selected
# - anomalies are computed relative to a mean over entire time series.
# - seasonal cycle is removed
# - a three point running mean is applied (the oceanic Nino index is 
#   defined as a 3-month running average of sst anomalies in the 
#   nino-3.4 region.)
# - a 12 point running mean function is also used
#
# - to use the data from MPI-ESM and HadGEM3, some regridding was necessary, using both
# - nco and cdo functions: 
# ncks -d time,0,1979 tos_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_18500116-20141216.nc test_1980ts.nc
# ncatted -a coordinates,tos,c,c,"latitude longitude" test_1980ts.nc
# cdo -L remapbil,mygrid -sethalo,-1,-1 test_1980ts.nc test_1980ts_latlon3.nc
#
# or:
# ncks -d time,0,1031 tos_Omon_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_20150116-21001216.nc test_ssp585ts.nc
#
# using this grid file: 
# gridtype = lonlat
# xsize    = 360
# ysize    = 40
# xfirst   = -179.5
# xinc     = 1
# yfirst   = -19.5
# yinc     = 1
#
# the ONI and RONI indices are computed using time series of the Nino3.4 and
# the tropical mean SSTs.  As shown in L'Heureux et al., 2024, the seasonal 
# cycle needs to be removed from both the Nino-3.4 and the tropical mean 
# becauase they have different seasonalities.  
#
# eventually we will want to calculate the mean of all the cmip timeseries. 
# will the zipper function be the best way to combine them before averaging?
# 
# functions:  fig_plot(), prepare_cmip_ts(), smooth_ts()
#
# cmip6 data was downloaded from the Climate Data Store:
# https://cds.climate.copernicus.eu/requests?tab=all
#
# levi silvers                                              sep 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using DataFrames
using CSV
using CairoMakie
using Statistics 
using NCDatasets

include("ensoFuncs.jl")

function fig_plot(inpv,xx,tit)
    #println(inpv)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(1000,200),
        )
    ax = Axis(f2[1,1];
        xlabel="time (monthly)",
        ylabel="temp anom (K)",
        title=tit,
        )
    #lines!(ax, 0.:1:xx-1, inpv -> cos(inpv);
    lines!(ax, 0.:1:xx-1, inpv;
    #scatterlines!(ax, xx, inpv;
        color=:black,
        linewidth=2,
        linestyle=:solid,
        )
    return f2
    #println(inpv[1:40])
end

#####

path="/Users/C823281551/"

# files with data from ssp585 are labelled with a b.   e.g.  file1b, file2b, file5b, file6b, and file9b
# CNRM
  # do we have ssp585 for the CNRM-CM6 model?  
#file1  = path*"data/cmip6/CNRMESM2/tos_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_18500116-20141216.nc" # --> doesn't appear to be regridded
file1  = path*"data/tos_CNRM_hist/tos_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_18500116-20141216.nc"
#file1b = path*"data/CNRM-CM6/tos_Omon_CNRM-CM6-1-HR_ssp585_r1i1p1f2_gn_20150116-21001216remapbil.nc"
#file1c = path*"data/cmip6/CNRMESM2/tos_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_20150116-21001216.nc" # --> doesn't appear to be regridded
#file1c = path*"data/CNRM-ESM2/tos_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_20150116-21001216.nc"
file1b = path*"data/cmip6/CNRMESM2/tos_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_20150116-21001216NewRegrid.nc"
# MPI-ESM
file2  = path*"data/MPI-ESM1/tos_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_18500116-20141216_regridded.nc"
file2b = path*"data/MPI-ESM1/tos_Omon_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_20150116-21001216_latlon.nc"
# GFDL
  # did GFDL not participate in the ssp simulations? 
file3  = path*"data/tos_GFDL_hist/tos_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_18500116-20141216.nc"
# E3SM
file4  = path*"data/E3SM/tos_Omon_E3SM-1-1-ECA_historical_r1i1p1f1_gr_18500116-20141216.nc"
file4b  = path*"data/"
# CESM2
#file5b = path*"data/CESM2/tos_Omon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215.nc" 
file5b = path*"data/CESM2/tos_Omon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215full_regrid.nc" 
# HadGEM3 (MB Andrews et al., 2020, JAMES?)
  # LL files have a nominal atmospheric resolution of 135km and an ocean resolution of 1 degrees.
  # MM files have a nominal atmospheric resolution of 60 km and an ocean resolution of 0.25 degrees.
file6  = path*"data/cmip6/HadGEM3/tos_Omon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_18500116-20141216_regridded.nc"
#file6b = path*"data/cmip6/HadGEM3/tos_Omon_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_20150116-21001216_regridded.nc"
file6b = path*"data/cmip6/HadGEM3/tos_Omon_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_20150116-21001216_regridFull.nc"
# ACCESS
file9b = path*"data/ACCESS-CM2/tos_Omon_ACCESS-CM2_ssp585_r1i1p1f1_gn_20150116-21001216_regridded.nc"
#FGOALS
file10b = path*"data/cmip6/FGOALS/tos_Omon_FGOALS-g3_ssp585_r1i1p1f1_gn_20150116-21001216_regrid.nc"

#--------------------------------
# incoming data in csv format:
file7 = path*"data/obs/observed_nino3.4.csv"
file8 = path*"data/obs/observed_tropicalmean.csv"
df1 = CSV.read(file7, header = 0, footerskip = 0, DataFrame) 
df2 = CSV.read(file8, header = 0, footerskip = 0, DataFrame) 
nms = ["year", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
dfa  = DataFrame(df1, nms)
dfb  = DataFrame(df2, nms)

  println("~~~~obs nino3.4~~~~~")
  println(collect(df1[1, 1:13]))
  println("~~~~obs trop mn~~~~~")
  println(collect(df2[1, 1:13]))

istart = 2
iend   = 170 # years
# this is necessary because of the weird structure of dataframes (1 year per row)
for i in istart:iend
    if i < istart + 1 
        global a1 = collect(df1[istart-1, 2:13]) # observed nino3.4
        global a2 = collect(df2[istart-1, 2:13]) # observed tr mn
    end
    b1 = collect(df1[i, 2:13]) # grab a row of the DataFrame and convert to vector
    b2 = collect(df2[i, 2:13]) # grab a row of the DataFrame and convert to vector
    global c1 = [a1; b1] # concatinate two vectors of observed nino3.4
    global c2 = [a2; b2] # concatinate two vectors of observed tropical mean
    a1 = c1
    a2 = c2
end

# compute the oceanic nino index (ONI) from nino3.4
#    nino34_ts_mn = mean(skipmissing(nino34_full))
enso3p4_mn = mean(c1)
nino3p4_anom = c1 .- enso3p4_mn    
tp_mn = mean(c2)
tp_anom = c2 .- tp_mn    

# remove the seasonal cycle
mn_oni = zeros(12)
#mn_oni = [mean(ts_oni[:,i]) for i in 2:13]
mn_a = [mean(df1[:,i]) for i in 2:13] # seasonal cycle of nino3.4
mn_b = [mean(df2[:,i]) for i in 2:13] # seasonal cycle of tropical mean
# one can check the seasonal cycle in the REPL with:
# lines(mn_a)
jend = 170*12

c1nsc  = zeros(jend)
c2nsc  = zeros(jend)
roni_a  = zeros(jend)

# remove the seasonal cycle before smoothing the time series.  
# ts_oni     -> nino3p4_anom
# tp_anom_sm -> tp_anom

for i in 1:12:jend
  # remove seasonal cycle of nino3p4
  c1nsc[i]=nino3p4_anom[i]-mn_a[1]
  c1nsc[i+1]=nino3p4_anom[i+1]-mn_a[2]
  c1nsc[i+2]=nino3p4_anom[i+2]-mn_a[3]
  c1nsc[i+3]=nino3p4_anom[i+3]-mn_a[4]
  c1nsc[i+4]=nino3p4_anom[i+4]-mn_a[5]
  c1nsc[i+5]=nino3p4_anom[i+5]-mn_a[6]
  c1nsc[i+6]=nino3p4_anom[i+6]-mn_a[7]
  c1nsc[i+7]=nino3p4_anom[i+7]-mn_a[8]
  c1nsc[i+8]=nino3p4_anom[i+8]-mn_a[9]
  c1nsc[i+9]=nino3p4_anom[i+9]-mn_a[10]
  c1nsc[i+10]=nino3p4_anom[i+10]-mn_a[11]
  c1nsc[i+11]=nino3p4_anom[i+11]-mn_a[12]
  # remove seasonal cycle of tropical mean 
  c2nsc[i]  =tp_anom[i]-mn_b[1]
  c2nsc[i+1]=tp_anom[i+1]-mn_b[2]
  c2nsc[i+2]=tp_anom[i+2]-mn_b[3]
  c2nsc[i+3]=tp_anom[i+3]-mn_b[4]
  c2nsc[i+4]=tp_anom[i+4]-mn_b[5]
  c2nsc[i+5]=tp_anom[i+5]-mn_b[6]
  c2nsc[i+6]=tp_anom[i+6]-mn_b[7]
  c2nsc[i+7]=tp_anom[i+7]-mn_b[8]
  c2nsc[i+8]=tp_anom[i+8]-mn_b[9]
  c2nsc[i+9]=tp_anom[i+9]-mn_b[10]
  c2nsc[i+10]=tp_anom[i+10]-mn_b[11]
  c2nsc[i+11]=tp_anom[i+11]-mn_b[12]
end

sig_oni   = std(c1nsc)   # standard deviation of oni
sig_dif   = std(c1nsc-c2nsc) # standard deviation of tr mean
sig_scale = sig_oni/sig_dif

# calculate a 3 point running mean
ts_oni = zeros(2040)
tmn_sm = zeros(2040)
istart= 2
jend  = 2040
for i in istart:jend-1
  ts_oni[i] = (c1nsc[i+1]+c1nsc[i]+c1nsc[i-1])/3
  tmn_sm[i] = (c2nsc[i+1]+c2nsc[i]+c2nsc[i-1])/3
end
ts_oni[1]    =ts_oni[2]
ts_oni[jend] =ts_oni[jend-1]
tmn_sm[1]    =tmn_sm[2]
tmn_sm[jend] =tmn_sm[jend-1]

roni_a = sig_scale.*(ts_oni - tmn_sm)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now deal with the data that is arriving from netcdf files: 
timelen = 1980
inpFile = file1
println("~~~~~~~~~~~~~~~~~~file 1~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
lat1 = 1
lat2 = 40
lon34a = 10
lon34b = 61
lat34a = 15
lat34b = 25

prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba1 = ts_rmn;
ba1nn = ts_rmn2;
###
inpFile = file2
println("~~~~~~~~~~~~~~~~~~file 2~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba2 = ts_rmn;
ba2nn = ts_rmn2;
##
inpFile = file3
println("~~~~~~~~~~~~~~~~~~file 3~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba3 = ts_rmn;
ba3nn = ts_rmn2;
#
inpFile = file4
println("~~~~~~~~~~~~~~~~~~file 4~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba4 = ts_rmn
#
inpFile = file6
println("~~~~~~~~~~~~~~~~~~file 6~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba6 = ts_rmn
#
##CNRM scenario timeseries...
timelen2=1032
inpFile = file1b # CNRM-ESM2
lat1 = 1
lat2 = 40
println("~~~~~~~~~~~~~~~~~~file 1b~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
lat1   = 71
lat2   = 110
lat34a = 85
lat34b = 96
lon34a = 10
lon34b = 61
prepare_cmip_ts(inpFile,timelen2,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba1b = ts_rmn
#prepare_cmip_ts(inpFile,timelen2,lat1,lat2,lon1,lon2)
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#ba1c = ts_rmn
#
##MPI scenario timeseries...
inpFile = file2b
println("~~~~~~~~~~~~~~~~~~file 2b~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
lat1 = 1
lat2 = 40
lon34a = 10
lon34b = 61
lat34a = 15
lat34b = 25
prepare_cmip_ts(inpFile,timelen2,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba2b = ts_rmn
#
inpFile = file5b # CESM2
println("~~~~~~~~~~~~~~~~~~file 5b~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
#prepare_cmip_ts(inpFile,timelen2,lat1,lat2,lon1,lon2)
lat1   = 71
lat2   = 110
lat34a = 85
lat34b = 96
lon34a = 10
lon34b = 61
prepare_cmip_ts(inpFile,timelen2,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba5b = ts_rmn
#
inpFile = file6b
println("~~~~~~~~~~~~~~~~~~file 6b~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)  # HadGEM3
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
lat1   = 71
lat2   = 110
lat34a = 85 
lat34b = 96 
lon34a = 10
lon34b = 61
prepare_cmip_ts(inpFile,timelen2,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba6b = ts_rmn
#
inpFile = file9b
println("~~~~~~~~~~~~~~~~~~file 9b~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile) # ACCESS
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
lat1   = 1
lat2   = 40
lat34a = 15 
lat34b = 26
lon34a = 10
lon34b = 61
prepare_cmip_ts(inpFile,timelen2,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba9b = ts_rmn
#
#FGOALS scenario timeseries...
timelen2=1032
inpFile = file10b # FGOALS
println("~~~~~~~~~~~~~~~~~~file 10b~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
lat1   = 71
lat2   = 110
lat34a = 85
lat34b = 96
lon34a = 10
lon34b = 61
prepare_cmip_ts(inpFile,timelen2,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba10b = ts_rmn
#prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
#
timelen = 1980
lat1 = 1
lat2 = 40
lon34a = 10
lon34b = 61
lat34a = 15
lat34b = 25
inpFile = file3
println("~~~~~~~~~~~~~~~~~~file 3~~~~~~~~~~~~~~~~~~~~~~")
println(inpFile)
println("""~~~~~~~~~~~~~~~~~~>>>>>>~~~~~~~~~~~~~~~~~~~~~~""")
prepare_cmip_ts(inpFile,timelen,lat1,lat2,lon34a,lon34b,lat34a,lat34b);
ba3 = ts_rmn
#
# define time axis for various experiment temporal ranges
A = collect(1854:1/12:2023.92)
B = collect(1850.0833333333333:1/12:2015)
C = collect(2015.083333:1/12:2101)
#D = collect(1990.083333:1/12:2023)

fig = Figure(;
    size = (800,300),
    )
ax = Axis(fig[1,1];
    xlabel="monthly mean, smoothed",
    ylabel="RONI anomalies",
    #xticks=([1850,1860,1870,1880,1890,1900,1910,1920,1930,1940,1950,1960,1970,1980,1990,2000,2010,2020]),
    xticks=([1850,1870,1890,1910,1930,1950,1970,1990,2010,2030,2050,2070,2090]),
    title="ENSO: historical period and the SSP5-8.5 scenario"
    )
smooth_12_ts(roni_a,2040)
blah3 = ts_12_sm
##lines!(ax, A,blah3[:], 
##    linewidth = 2.,
##    color = "black",
##    label = "Observed: RONI"
##    )
##limits!(1850, 2100, -4, 4)
#
## historical
#
smooth_12_ts(ba1,timelen)
ba1_sm = ts_12_sm
##smooth_12_ts(ba1nn,timelen)
##ba1nn_sm = ts_12_sm
lines!(ax, B,ba1_sm[:], 
    linewidth = 0.75,
    color = "grey77"
    #color = "paleturquoise1"
    #label = "CNRM: RONI"
    )
limits!(1850, 2100, -4, 4)
#
## ssp585
#
smooth_12_ts(ba1b,timelen2)
ba1b_sm = ts_12_sm
lines!(ax, C,ba1b_sm[:], 
    linewidth = 0.75,
    color = "grey77"
    #color = "paleturquoise2"
    #label = CNRMESM2
    )
#smooth_12_ts(ba1c,timelen2)
#ba1c_sm = ts_12_sm
#lines!(ax, C,ba1c_sm[:], 
#    linewidth = 0.75,
#    color = "paleturquoise3"
#    #label = "CNRM: RONI"
#    )
#
smooth_12_ts(ba2,timelen)
ba2_sm = ts_12_sm #  historical
smooth_12_ts(ba2b,timelen2)
ba2b_sm = ts_12_sm # ssp585
##smooth_12_ts(ba2nn,timelen)
##ba2nn_sm = ts_12_sm
lines!(ax, B,ba2_sm[:], 
    linewidth = 0.75,
    color = "grey77"
    #color = "aquamarine"
    #label = "MPI: RONI"
    )
lines!(ax, C,ba2b_sm[:], 
    linewidth = 0.75,
    color = "grey77"
    #color = "aquamarine"
    #label = " MPI-ESM1: RONI"
    )
smooth_12_ts(ba3,timelen)
ba3_sm = ts_12_sm # historical
##smooth_12_ts(ba3nn,timelen)
##ba3nn_sm = ts_12_sm
lines!(ax, B,ba3_sm[:], 
    linewidth = 0.75,
    color = "grey77"
    #color = "lightcyan"
    #label = "GFDL: RONI"
    )

smooth_12_ts(ba4,timelen)
ba4_sm = ts_12_sm # historical
lines!(ax, B,ba4_sm[:], 
    linewidth = 0.75,
    color = "grey77"
    #color = "paleturquoise"
    #label = "E3SM: RONI"
    )
#
smooth_12_ts(ba5b,timelen2)
ba5b_sm = ts_12_sm # ssp585
lines!(ax, C,ba5b_sm[:], 
    linewidth = 0.75,
    color = "grey77"
    #color = "mistyrose"
    #label = "CESM2: RONI"
    )
#
smooth_12_ts(ba6,timelen)
ba6_sm = ts_12_sm # historical
smooth_12_ts(ba6b,timelen2)
ba6b_sm = ts_12_sm # ssp585
lines!(ax, B,ba6_sm[:], 
    linewidth = 0.75,
    color = "grey77"
    #color = "thistle2"
    #label = "HadGEM3: RONI"
    )
lines!(ax, C,ba6b_sm[:], 
    linewidth = 0.75,
    color = "grey77"
    #color = "thistle"
    #label = "HadGEM3: RONI"
    )
#
smooth_12_ts(ba9b,timelen2)
ba9b_sm = ts_12_sm # ssp585
lines!(ax, C,ba9b_sm[:], 
    linewidth = 0.75,
    color = "grey77"
    #color = "lavender"
    #label = "ACCESS: RONI"
    )
#
smooth_12_ts(ba10b,timelen2)
ba10b_sm = ts_12_sm # ssp585
lines!(ax, C,ba10b_sm[:], 
    linewidth = 0.75,
    color = "grey77"
    #color = "black"
    #label = "ACCESS: RONI"
    )
#
## ssp585
tenCent  = [ba1b_sm';ba2b_sm';ba5b_sm';ba9b_sm';ba6b_sm']
mnSSP = mean(tenCent, dims = 1) # i think that averaging these timeseries results in a matrix
# which then is difficult to use in coordination with an array.   

# plot individual model for ssp585 if desired.
#mnSSP = ba6b_sm

## historical
lightOut = [ba1_sm';ba2_sm';ba3_sm';ba4_sm';ba6_sm']
##lightOut = [ba2_sm';ba3_sm';ba4_sm';ba6_sm']
mnH   = mean(lightOut, dims = 1)
#
# compute the trends
numyrs = 33; # number of years of observations to compute the trend over. 
firstyr=1990.;
lastyr=12*numyrs-1;
#lastyr=12*33-1; # this was the default

i1=firstyr+0.083333
i2=firstyr+.083333333333333+lastyr/12

#timeScat = collect(0:12*33-1);
timeScat = collect(0:lastyr);
#D = collect(1990.083333:1/12:2023)
D = collect(i1:1/12:i2)
lD = size(D)
lengthD = lD[1]
# for SSP5-8.5 simulations
timeB = collect(0:1031.);

# roni_a --> 2040 elements, unsmoothed roni timeseries
#            the timeseries extends from Jan 1854 to Dec 2023

# works! histTrend = find_best_fit(timeA,ba1_sm);
#histTrend = find_best_fit(timeScat,mnH[1585:1585+lastyr]);

#histTrend = find_best_fit(timeScat,roni_a[1585:1585+lastyr]);
histTrend = find_best_fit(timeScat,roni_a[end-lengthD+1:end]);

# below computes trend of historical model simulation
#histTrend = find_best_fit(timeScat,mnH[1585:1980]);
# project the trend to the correct timeslice using the least squares coefficients:
hTrend    = collect(timeScat).*histTrend[1] .+ histTrend[2];

#cmipTrend = find_best_fit(timeB,mnSSP);
cmipTrend = find_best_fit(timeB,mnSSP[:]);
cmTrend   = collect(timeB).*cmipTrend[1] .+ cmipTrend[2];

##smooth_12_ts(ba6,timelen)
##ba6_sm = ts_12_sm # historical
lines!(ax, B,mnH[:], 
    linewidth = 2.0,
    color = "red",
    label = "mn CMIP6"
    )
#lines!(ax, B,hTrend)
lines!(ax, C,mnSSP[:], 
    linewidth = 2.0,
    color = "red"
    #label = "mn SSP"
    )
lines!(ax, C,cmTrend,
    color = "red",
    linewidth = 3.0,
    )
lines!(ax, A,blah3[:], 
#lines!(ax, A,roni_a[:], 
    linewidth = 2.,
    color = "black",
    label = "Observed"
    )
lines!(ax, D,hTrend,
    linewidth = 3.0,
    color = "black"
    )
#
axislegend( position=:lt)
#
save("plotENSOcmip_testA.png",fig)








