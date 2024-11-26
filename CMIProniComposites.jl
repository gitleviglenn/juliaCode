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
# both SST (tos) and RH (hur) fields need to be read into script
#
# this should be useful for creating plots with EArth data: 
# https://forem.julialang.org/gage/creating-nice-simple-geo-plots-from-scratch-using-makiejl-301
# 
# levi silvers                                                                nov 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using DataFrames
using CairoMakie
using Statistics 
using NCDatasets

# calculate RONI
# based on a threshold of RONI determine ENSO events
# grab RH fields when ENSO events are occuring
# plot the composite RH fields

function check_thresh_high(inpFile, inpTS, thshd)
    ds = NCDataset(inpFile)
    ds.attrib
    nctime = ds["time"]
    println("~~~~~~~~~Silence~~~~~~~~~~~~~~~~~~")
    println(nctime[1:10])
    println("~~~~~~~~~Golden~~~~~~~~~~~~~~~~~~")
    # for incoming file, get the time dimension 
    # for the incoming timeseries, check at each index, if the value is greater than 
    # the treshold, if so return true, otherwise return false.   
    # or perhaps return 1 or 0.
    # for i in full lenght, is val >= thshd, grab i value
    #
    global high = findall(inpTS .> thshd)  
    #
    println("~~~~~~~~~high Threshold values~~~~~~~~~~~~~~~~~~")
    println(inpTS[high])
    println("~~~~~~~~~End high Threshold Values~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~high Threshold Times     ~~~~~~~~~~~~~~~~~~")
    println(nctime[high])
    println("~~~~~~~~~End high Threshold Times     ~~~~~~~~~~~~~~~~~~")
end

function check_thresh_low(inpFile, inpTS, thshd)
    ds = NCDataset(inpFile)
    ds.attrib
    nctime = ds["time"]
    # for incoming file, get the time dimension 
    # for the incoming timeseries, check at each index, if the value is greater than 
    # the treshold, if so return true, otherwise return false.   
    # or perhaps return 1 or 0.
    # for i in full lenght, is val >= thshd, grab i value
    #
    global low = findall(inpTS .< thshd)  
    #
    println("~~~~~~~~~low Threshold values~~~~~~~~~~~~~~~~~~")
    println(inpTS[low])
    println("~~~~~~~~~End low Threshold Values~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~low Threshold Times     ~~~~~~~~~~~~~~~~~~")
    println(nctime[low])
    println("~~~~~~~~~low End Threshold Times     ~~~~~~~~~~~~~~~~~~")
end

function prepare_cmip_ts(inpFile,len)
  #
  # checks for missing values and NaNs.
  # removes the seasonal cycle
  # compute a 3 point running mean
  #
    ds1 = NCDataset(inpFile)
    ds1.attrib
    # to print meta data for a particular variable: 
    ds1["tos"]
    sst1 = ds1["tos"]
    nclat = ds1["lat"]
    nclon = ds1["lon"]
    nctime = ds1["time"]
    println(nclat[15:26])
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~long points: ~~~~~~~~~~~~~~~~")
    println(nclon[10:61])
    println("~~~~~lat points: ~~~~~~~~~~~~~~~~")
    println(nclat[15:26])
    #
    nino34_full = sst1[10:61, 15:26, :];
    trop_full   = sst1[:, :, :];
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~type of nino34_full ~~~~~~~~~~~~~~~~~~~~~~~")
    typeof(nino34_full)
    size(nino34_full)
    println(nino34_full[:,7,10])
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    
    #tlength = 1980
    tlength = len 
    nino34_ts = zeros(tlength)
    tropmn_ts = zeros(tlength)
    
    nino34_ts_mn = mean(filter(!isnan, skipmissing(nino34_full)))
    tropmn_ts_mn = mean(filter(!isnan, skipmissing(trop_full)))
    
    #println("fury fury fury")
    #println(tropmn_ts_mn) 
    # compute the anomalies
    for i = 1:tlength
        #nino34_ts[i]=mean(skipmissing(nino34_full[:,:,i]))-nino34_ts_mn
        nino34_ts[i]=mean(filter(!isnan, skipmissing(nino34_full[:,:,i])))-nino34_ts_mn
        tropmn_ts[i]=mean(filter(!isnan, skipmissing(trop_full[:,:,i])))-tropmn_ts_mn
    end
    
    # remove the seasonal cycle:
    # note: these months are not necessarily correct.  it depends on what month the ts starts on...
    ss  = zeros(12)
    jan = [mean([nino34_ts[i] for i in 1:12:tlength])]
    feb = [mean([nino34_ts[i] for i in 2:12:tlength])]
    mar = [mean([nino34_ts[i] for i in 3:12:tlength])]
    apr = [mean([nino34_ts[i] for i in 4:12:tlength])]
    may = [mean([nino34_ts[i] for i in 5:12:tlength])]
    jun = [mean([nino34_ts[i] for i in 6:12:tlength])]
    jul = [mean([nino34_ts[i] for i in 7:12:tlength])]
    aug = [mean([nino34_ts[i] for i in 8:12:tlength])]
    sep = [mean([nino34_ts[i] for i in 9:12:tlength])]
    oct = [mean([nino34_ts[i] for i in 10:12:tlength])]
    nov = [mean([nino34_ts[i] for i in 11:12:tlength])]
    dec = [mean([nino34_ts[i] for i in 12:12:tlength])]
    
    ss = [jan feb mar apr may jun jul aug sep oct nov dec]

    global ts_rmn_nsc = zeros(tlength)
    for i in 1:12:tlength
      ts_rmn_nsc[i]    = nino34_ts[i]    - ss[1]
      ts_rmn_nsc[i+1]  = nino34_ts[i+1]  - ss[2]
      ts_rmn_nsc[i+2]  = nino34_ts[i+2]  - ss[3]
      ts_rmn_nsc[i+3]  = nino34_ts[i+3]  - ss[4]
      ts_rmn_nsc[i+4]  = nino34_ts[i+4]  - ss[5]
      ts_rmn_nsc[i+5]  = nino34_ts[i+5]  - ss[6]
      ts_rmn_nsc[i+6]  = nino34_ts[i+6]  - ss[7]
      ts_rmn_nsc[i+7]  = nino34_ts[i+7]  - ss[8]
      ts_rmn_nsc[i+8]  = nino34_ts[i+8]  - ss[9]
      ts_rmn_nsc[i+9]  = nino34_ts[i+9]  - ss[10]
      ts_rmn_nsc[i+10] = nino34_ts[i+10] - ss[11]
      ts_rmn_nsc[i+11] = nino34_ts[i+11] - ss[12]
    end

    ss1  = zeros(12)
    jan = [mean([tropmn_ts[i] for i in 1:12:tlength])]
    feb = [mean([tropmn_ts[i] for i in 2:12:tlength])]
    mar = [mean([tropmn_ts[i] for i in 3:12:tlength])]
    apr = [mean([tropmn_ts[i] for i in 4:12:tlength])]
    may = [mean([tropmn_ts[i] for i in 5:12:tlength])]
    jun = [mean([tropmn_ts[i] for i in 6:12:tlength])]
    jul = [mean([tropmn_ts[i] for i in 7:12:tlength])]
    aug = [mean([tropmn_ts[i] for i in 8:12:tlength])]
    sep = [mean([tropmn_ts[i] for i in 9:12:tlength])]
    oct = [mean([tropmn_ts[i] for i in 10:12:tlength])]
    nov = [mean([tropmn_ts[i] for i in 11:12:tlength])]
    dec = [mean([tropmn_ts[i] for i in 12:12:tlength])]
    
    ss1 = [jan feb mar apr may jun jul aug sep oct nov dec]

    global ts_trmn_nsc = zeros(tlength)
    for i in 1:12:tlength
      ts_trmn_nsc[i]    = tropmn_ts[i]    - ss1[1]
      ts_trmn_nsc[i+1]  = tropmn_ts[i+1]  - ss1[2]
      ts_trmn_nsc[i+2]  = tropmn_ts[i+2]  - ss1[3]
      ts_trmn_nsc[i+3]  = tropmn_ts[i+3]  - ss1[4]
      ts_trmn_nsc[i+4]  = tropmn_ts[i+4]  - ss1[5]
      ts_trmn_nsc[i+5]  = tropmn_ts[i+5]  - ss1[6]
      ts_trmn_nsc[i+6]  = tropmn_ts[i+6]  - ss1[7]
      ts_trmn_nsc[i+7]  = tropmn_ts[i+7]  - ss1[8]
      ts_trmn_nsc[i+8]  = tropmn_ts[i+8]  - ss1[9]
      ts_trmn_nsc[i+9]  = tropmn_ts[i+9]  - ss1[10]
      ts_trmn_nsc[i+10] = tropmn_ts[i+10] - ss1[11]
      ts_trmn_nsc[i+11] = tropmn_ts[i+11] - ss1[12]
    end

    # smooth the time series with a running mean
    global ts_nino34_rmn = zeros(tlength)
    global ts_tropmn_rmn = zeros(tlength)
    global ts_rmn        = zeros(tlength)
    global ts_rmn2        = zeros(tlength)
    istart= 2
    jend  = tlength
    for i in istart:jend-1
      ts_nino34_rmn[i] = (ts_rmn_nsc[i+1]+ts_rmn_nsc[i]+ts_rmn_nsc[i-1])/3
      ts_tropmn_rmn[i] = (ts_trmn_nsc[i+1]+ts_trmn_nsc[i]+ts_trmn_nsc[i-1])/3
    end
    ts_nino34_rmn[1]=ts_nino34_rmn[2]
    ts_nino34_rmn[jend]=ts_nino34_rmn[jend-1]
    ts_tropmn_rmn[1]=ts_tropmn_rmn[2]
    ts_tropmn_rmn[jend]=ts_tropmn_rmn[jend-1]

    # compute the scaling factor
    sig_oni   =  std(ts_nino34_rmn)
    sig_diff  =  std(ts_nino34_rmn-ts_tropmn_rmn)
    sig_scale =  sig_oni/sig_diff

    #print(ts_nino34_rmn)
    #println("fury and hatred")
    # compute the relative oni index (RONI):
    ts_rmn = sig_scale.*(ts_nino34_rmn-ts_tropmn_rmn)
    ts_rmn2 = sig_scale.*(ts_nino34_rmn)
    return ts_rmn 
end


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# beginning of work with GCM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path="/Users/C823281551/data/"

# historical
#file1  = path*
modelp="MPI-ESM1-2-LR"
# SSP585 scenario
filehur  = path*"cmip6/CESM2/hur_Amon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215.nc"
file1b  = path*"cmip6/CESM2/tos_Omon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215.nc" 
#filehur  = path*"cmip6/MPIESM/hur_Amon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216.nc"
#file1b  = path*"cmip6/MPIESM/tos_Omon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_latlon.nc" 

data  = NCDataset(filehur)

data.attrib

rh   = data["hur"] # hur(time, plev, lat, lon)
lat  = data["lat"]
lon  = data["lon"]
tme = data["time"]
lev  = data["plev"]

tim = 5 # we need to grab rh values at the times selected by check_thresh
rh1 = rh[:,:,1,tim]
rh2 = rh[:,:,2,tim]


#scenario timeseries
timelen2=1032
inpFile = file1b
prepare_cmip_ts(inpFile,timelen2)
ba1b = ts_rmn

# smooth time series to be checked for threshold
smooth_12_ts(ba1b,timelen2)
ba1b_sm1 = ts_12_sm

thshd = 2.0
inpFile = file1b
# i think the only reason we need to pass in a file is to get the 
# time dimension from it.   to speed things up we could pass that
# in instead of dealing with the file within check_thresh()
check_thresh_high(inpFile, ba1b_sm1, thshd) # output is 'high'
thshd = -1.5
check_thresh_low(inpFile, ba1b_sm1, thshd)  # output is 'low'

# we need output times from check_thresh() at which to grab RH vals.
#
# read in list of times.   for an index, i, check if time[i] is in
# the incoming list.   If so, grap RH[i]

ds = NCDataset(filehur)
ds.attrib
nctime = ds["time"]
println("~~~~~~~~~Silence~~~~~~~~~~~~~~~~~~")
println(nctime[high[1:10]])
println("~~~~~~~~~Golden~~~~~~~~~~~~~~~~~~~")

#rhhigh = rh[:,:,1,high]

#test_rh = Array{Float64, 3}(undef, 288, 84, 2)

## instead of 2 the last dimension should be the number of values in high (size(high))
#test_rh = Array{Union{Missing, Float64}, 3}(undef, 288, 84, 2)
#rh_high = Array{Union{Missing, Float64}, 3}(undef, 288, 84, 10)
#rh_low  = Array{Union{Missing, Float64}, 3}(undef, 288, 84, 10)
dims = size(rh)
test_rh = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 2)
#test_rh = Array{Union{Missing, Float64}, 3}(undef, 192, 96, 2)
rh_high = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 10)
rh_low  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 10)

#test_rh = Matrix{Float64}(undef, 288, 84)
test_rh[:,:,1] = rh[:,:,1,high[1]]
#test_rh2 = Matrix{Float64}(undef, 288, 84)
test_rh[:,:,2] = rh[:,:,1,high[2]]
# restults can be plotted in the REPL as: 
# contourf(test_rh[:,:,1])

rh_high[:,:,1] = rh[:,:,1,high[1]]
rh_high[:,:,2] = rh[:,:,1,high[2]]
rh_high[:,:,3] = rh[:,:,1,high[3]]
rh_high[:,:,4] = rh[:,:,1,high[4]]
rh_high[:,:,5] = rh[:,:,1,high[5]]
rh_high[:,:,6] = rh[:,:,1,high[6]]
rh_high[:,:,7] = rh[:,:,1,high[7]]
rh_high[:,:,8] = rh[:,:,1,high[8]]
rh_high[:,:,9] = rh[:,:,1,high[9]]
rh_high[:,:,10] = rh[:,:,1,high[10]]

rh_low[:,:,1]  = rh[:,:,1,low[1]]
rh_low[:,:,2]  = rh[:,:,1,low[2]]
rh_low[:,:,3]  = rh[:,:,1,low[3]]
rh_low[:,:,4]  = rh[:,:,1,low[4]]
rh_low[:,:,5]  = rh[:,:,1,low[5]]
rh_low[:,:,6]  = rh[:,:,1,low[6]]
rh_low[:,:,7]  = rh[:,:,1,low[7]]
rh_low[:,:,8]  = rh[:,:,1,low[8]]
rh_low[:,:,9]  = rh[:,:,1,low[9]]
rh_low[:,:,10]  = rh[:,:,1,low[10]]

rh_high_mn = mean(rh_high, dims = 3)
rh_low_mn  = mean(rh_low, dims = 3)

rh_diff = rh_high_mn - rh_low_mn

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
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
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
function fig_anom_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(600,300),
        )
    ax = Axis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        xlabel="longitude",
        ylabel="latitude",
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             levels = range(-20, 20, length = 100), 
             #colormap = :batlow,
             colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        Colorbar(f2[1,2], bb)
    return f2
end

rh_2_plot = rh_low_mn[:,:,1]
fig = fig_1_plot(rh_2_plot,lon,lat,"RH low")
#rh_2_plot = rh_high_mn[:,:,1]
#fig = fig_1_plot(rh_2_plot,lon,lat,"RH high")

rh_anom = rh_diff[:,:,1]
fig = fig_anom_plot(rh_anom,lon,lat,"RH anomaly, MPI ssp585")


