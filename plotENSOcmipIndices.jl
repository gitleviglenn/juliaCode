#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotENSOcmipIndices.jl
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

function prepare_cmip_ts(inpFile,len,troplat1,troplat2,ln1,ln2,lt1,lt2)
#function prepare_cmip_ts(inpFile,len)
    # the seasonal cycle needs to be removed
    # the tropical mean in L'Heureux are defined as +/-20, it does not 
    # appear that is being taken into account for potentially different grids.  Check!
    ds1 = NCDataset(inpFile)
    ds1.attrib
    # to print meta data for a particular variable: 
    ds1["tos"]
    sst1 = ds1["tos"] # what are the dimensions here?   Should be +/-20 degrees
    nclat = ds1["lat"]
    nclon = ds1["lon"]
    nctime = ds1["time"]
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~ nino34 long points: ~~~~~~~~~~~~~~~~")
    #println(nclon[10:61])
    println(nclon[ln1:ln2])
    println("~~~~~ nino34 lat points: ~~~~~~~~~~~~~~~~")
    println(nclat[lt1:lt2])
    println("~~~~~tropical boundary lat points: ~~~~~~~~~~~~~~~~")
    println(nclat[troplat1:troplat2])
    #
    nino34_full = sst1[ln1:ln2, lt1:lt2, :];
    trop_mean   = sst1[:, troplat1:troplat2, :];
    #trop_mean   = sst1[:, :, :];
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    #println("~~~~type of nino34_full ~~~~~~~~~~~~~~~~~~~~~~~")
    #typeof(nino34_full)
    #size(nino34_full)
    #println(nino34_full[:,7,10])
    #println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    
    tlength = len 
    global nino34_ts = zeros(tlength)
    tropmn_ts = zeros(tlength)
   
    # compute climatological average value in nino34 and tropical regions 
    nino34_ts_mn = mean(filter(!isnan, skipmissing(nino34_full)))
    tropmn_ts_mn = mean(filter(!isnan, skipmissing(trop_mean)))
    
    # compute the anomalies
    # nino3.4 index --> nino34_ts
    # ONI           --> 3 month runnign mean of nino34_ts
    for i = 1:tlength
        #nino34_ts[i]=mean(skipmissing(nino34_full[:,:,i]))-nino34_ts_mn
        nino34_ts[i]=mean(filter(!isnan, skipmissing(nino34_full[:,:,i])))-nino34_ts_mn;
        tropmn_ts[i]=mean(filter(!isnan, skipmissing(trop_mean[:,:,i])))-tropmn_ts_mn;
    end
    
    # remove the seasonal cycle:
    # note: these months are not necessarily correct.  it depends on what month the ts starts on...
    ss    = zeros(12)
    sst   = zeros(12)
    jan   = [mean([nino34_ts[i] for i in 1:12:tlength])]
    janst = std(nino34_ts[i] for i in 1:12:tlength)
    feb   = [mean([nino34_ts[i] for i in 2:12:tlength])]
    febst = std(nino34_ts[i] for i in 2:12:tlength)
    mar   = [mean([nino34_ts[i] for i in 3:12:tlength])]
    marst = std(nino34_ts[i] for i in 3:12:tlength)
    apr   = [mean([nino34_ts[i] for i in 4:12:tlength])]
    aprst = std(nino34_ts[i] for i in 4:12:tlength)
    may   = [mean([nino34_ts[i] for i in 5:12:tlength])]
    mayst = std(nino34_ts[i] for i in 5:12:tlength)
    jun   = [mean([nino34_ts[i] for i in 6:12:tlength])]
    junst = std(nino34_ts[i] for i in 6:12:tlength)
    jul   = [mean([nino34_ts[i] for i in 7:12:tlength])]
    julst = std(nino34_ts[i] for i in 7:12:tlength)
    aug   = [mean([nino34_ts[i] for i in 8:12:tlength])]
    augst = std(nino34_ts[i] for i in 8:12:tlength)
    sep   = [mean([nino34_ts[i] for i in 9:12:tlength])]
    sepst = std(nino34_ts[i] for i in 9:12:tlength)
    oct   = [mean([nino34_ts[i] for i in 10:12:tlength])]
    octst = std(nino34_ts[i] for i in 10:12:tlength)
    nov   = [mean([nino34_ts[i] for i in 11:12:tlength])]
    novst = std(nino34_ts[i] for i in 11:12:tlength)
    dec   = [mean([nino34_ts[i] for i in 12:12:tlength])]
    decst = std(nino34_ts[i] for i in 12:12:tlength)
    
    ss    = [jan feb mar apr may jun jul aug sep oct nov dec]
    sst   = [janst febst marst aprst mayst junst julst augst sepst octst novst decst]

    # computed ts for nino 3.4 minus seasonal cycle
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

#-------------------------------------------------------------------------- 
    # calculate ONI: compute the 3mn running average
    global ts_oni_cmip        = zeros(tlength)
    istart= 2
    jend  = tlength
    for i in istart:jend-1
      ts_oni_cmip[i] = (ts_rmn_nsc[i+1]+ts_rmn_nsc[i]+ts_rmn_nsc[i-1])/3;
    end
    ts_oni_cmip[1]=ts_oni_cmip[2]
    ts_oni_cmip[jend]=ts_oni_cmip[jend-1]

    # now compute the standard deviation of ONI
    st_oni    = zeros(12)
    janst_oni = std(ts_oni_cmip[i] for i in 1:12:tlength)
    febst_oni = std(ts_oni_cmip[i] for i in 2:12:tlength)
    marst_oni = std(ts_oni_cmip[i] for i in 3:12:tlength)
    aprst_oni = std(ts_oni_cmip[i] for i in 4:12:tlength)
    mayst_oni = std(ts_oni_cmip[i] for i in 5:12:tlength)
    junst_oni = std(ts_oni_cmip[i] for i in 6:12:tlength)
    julst_oni = std(ts_oni_cmip[i] for i in 7:12:tlength)
    augst_oni = std(ts_oni_cmip[i] for i in 8:12:tlength)
    sepst_oni = std(ts_oni_cmip[i] for i in 9:12:tlength)
    octst_oni = std(ts_oni_cmip[i] for i in 10:12:tlength)
    novst_oni = std(ts_oni_cmip[i] for i in 11:12:tlength)
    decst_oni = std(ts_oni_cmip[i] for i in 12:12:tlength)
    st_oni    = [janst_oni febst_oni marst_oni aprst_oni mayst_oni junst_oni julst_oni augst_oni sepst_oni octst_oni novst_oni decst_oni]

    global ts_oni_st = zeros(tlength)
    for i in 1:12:tlength
      ts_oni_st[i]     = st_oni[1]
      ts_oni_st[i+1]   = st_oni[2]
      ts_oni_st[i+2]   = st_oni[3]
      ts_oni_st[i+3]   = st_oni[4]
      ts_oni_st[i+4]   = st_oni[5]
      ts_oni_st[i+5]   = st_oni[6]
      ts_oni_st[i+6]   = st_oni[7]
      ts_oni_st[i+7]   = st_oni[8]
      ts_oni_st[i+8]   = st_oni[9]
      ts_oni_st[i+9]   = st_oni[10]
      ts_oni_st[i+10]  = st_oni[11]
      ts_oni_st[i+11]  = st_oni[12]
    end
    # ts_oni_st is the standard deviation of ONI: sigma(oni)
#-------------------------------------------------------------------------- 


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

    # compute ts for the tropical mean sst minus seasonal cycle
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
    global ts_rmn2       = zeros(tlength)
    global oniMtpm       = zeros(tlength)
    istart= 2
    jend  = tlength
    for i in istart:jend-1
      ts_nino34_rmn[i] = (ts_rmn_nsc[i+1]+ts_rmn_nsc[i]+ts_rmn_nsc[i-1])/3;
      ts_tropmn_rmn[i] = (ts_trmn_nsc[i+1]+ts_trmn_nsc[i]+ts_trmn_nsc[i-1])/3;
    end
    ts_nino34_rmn[1]=ts_nino34_rmn[2]
    ts_nino34_rmn[jend]=ts_nino34_rmn[jend-1]
    ts_tropmn_rmn[1]=ts_tropmn_rmn[2]
    ts_tropmn_rmn[jend]=ts_tropmn_rmn[jend-1]

    #--------------------------------------------------
    # compute the unscaled RONI
    oniMtpm = ts_oni_cmip .- ts_tropmn_rmn
    #--------------------------------------------------

    #--------------------------------------------------
    # compute the standard deviation of oni - tropmn
    st_oniMtm    = zeros(12)
    janst_oniMtm = std(oniMtpm[i] for i in 1:12:tlength)
    febst_oniMtm = std(oniMtpm[i] for i in 2:12:tlength)
    marst_oniMtm = std(oniMtpm[i] for i in 3:12:tlength)
    aprst_oniMtm = std(oniMtpm[i] for i in 4:12:tlength)
    mayst_oniMtm = std(oniMtpm[i] for i in 5:12:tlength)
    junst_oniMtm = std(oniMtpm[i] for i in 6:12:tlength)
    julst_oniMtm = std(oniMtpm[i] for i in 7:12:tlength)
    augst_oniMtm = std(oniMtpm[i] for i in 8:12:tlength)
    sepst_oniMtm = std(oniMtpm[i] for i in 9:12:tlength)
    octst_oniMtm = std(oniMtpm[i] for i in 10:12:tlength)
    novst_oniMtm = std(oniMtpm[i] for i in 11:12:tlength)
    decst_oniMtm = std(oniMtpm[i] for i in 12:12:tlength)
    st_oniMtm    = [janst_oniMtm febst_oniMtm marst_oniMtm aprst_oniMtm mayst_oniMtm junst_oniMtm julst_oniMtm augst_oniMtm sepst_oniMtm octst_oniMtm novst_oniMtm decst_oniMtm]

    global tsoniMtm_st = zeros(tlength)
    for i in 1:12:tlength
      tsoniMtm_st[i]     = st_oniMtm[1]
      tsoniMtm_st[i+1]   = st_oniMtm[2]
      tsoniMtm_st[i+2]   = st_oniMtm[3]
      tsoniMtm_st[i+3]   = st_oniMtm[4]
      tsoniMtm_st[i+4]   = st_oniMtm[5]
      tsoniMtm_st[i+5]   = st_oniMtm[6]
      tsoniMtm_st[i+6]   = st_oniMtm[7]
      tsoniMtm_st[i+7]   = st_oniMtm[8]
      tsoniMtm_st[i+8]   = st_oniMtm[9]
      tsoniMtm_st[i+9]   = st_oniMtm[10]
      tsoniMtm_st[i+10]  = st_oniMtm[11]
      tsoniMtm_st[i+11]  = st_oniMtm[12]
    end
    # tsoniMtm_st is the standard deviation of ONI-tropmn: sigma(oni-tmn)
    #--------------------------------------------------

    # compute the scaling factor
    sig_oni   =  std(ts_nino34_rmn)
    sig_diff  =  std(ts_nino34_rmn-ts_tropmn_rmn)
    sig_scale =  sig_oni/sig_diff

    sigma     =  ts_oni_st./tsoniMtm_st

    #print(ts_nino34_rmn)
    #println("fury and hatred")
    # compute the relative oni index (RONI):
    #ts_rmn = sig_scale.*(ts_nino34_rmn-ts_tropmn_rmn);
    ts_rmn = sigma.*oniMtpm;
    ts_rmn2 = sig_scale.*(ts_nino34_rmn);
    return ts_rmn;
end

function smooth_ts(inpTS,len)
    # smooth the time series with a running mean
    global ts_sm = zeros(len)
    istart= 2
    jend  = len
    for i in istart:jend-1
      ts_sm[i] = (inpTS[i+1]+inpTS[i]+inpTS[i-1])/3
    end
    ts_sm[1]=ts_sm[2]
    ts_sm[jend]=ts_sm[jend-1]

    return ts_sm
end

function smooth_12_ts(inpTS,len)
    # smooth the time series with a running mean
    global ts_12_sm = zeros(len)
    istart= 6
    jend  = len
    for i in istart:jend-6
      ts_12_sm[i] = (inpTS[i+6]+inpTS[i+5]+inpTS[i+4]+inpTS[i+3]+inpTS[i+2]+inpTS[i+1]+inpTS[i]+inpTS[i-1]+inpTS[i-2]+inpTS[i-3]+inpTS[i-4]+inpTS[i-5])/12
    end
    return ts_12_sm
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
iend   = 170
for i in istart:iend
    if i < istart + 1 
        global a1 = collect(df1[istart-1, 2:13]) # observed nino3.4
        global a2 = collect(df2[istart-1, 2:13]) # observed tr mn
    end
    b1 = collect(df1[i, 2:13]) # grab a row of the DataFrame and convert to vector
    b2 = collect(df2[i, 2:13]) # grab a row of the DataFrame and convert to vector
    global c1 = [a1; b1] # concatinate two vectors
    global c2 = [a2; b2] # concatinate two vectors
    a1 = c1
    a2 = c2
end

# compute the oceanic nino index from nino3.4
#    nino34_ts_mn = mean(skipmissing(nino34_full))
enso3p4_mn = mean(c1)
nino3p4_anom = c1 .- enso3p4_mn    
tp_mn = mean(c2)
tp_anom = c2 .- tp_mn    

# remove the seasonal cycle
mn_oni = zeros(12)
#mn_oni = [mean(ts_oni[:,i]) for i in 2:13]
mn_a = [mean(df1[:,i]) for i in 2:13]
mn_b = [mean(df2[:,i]) for i in 2:13]
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
#
timelen = 1980
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

fig = Figure(;
    size = (800,300),
    )
ax = Axis(fig[1,1];
    xlabel="monthly mean, smoothed",
    ylabel="ENSO RONI anomalies",
    #xticks=([1850,1860,1870,1880,1890,1900,1910,1920,1930,1940,1950,1960,1970,1980,1990,2000,2010,2020]),
    xticks=([1850,1870,1890,1910,1930,1950,1970,1990,2010,2030,2050,2070,2090]),
    title="ENSO: historical period and ssp585"
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
    color = "paleturquoise1"
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
    color = "paleturquoise2"
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
    color = "aquamarine"
    #label = "MPI: RONI"
    )
lines!(ax, C,ba2b_sm[:], 
    linewidth = 0.75,
    color = "aquamarine"
    #label = " MPI-ESM1: RONI"
    )
smooth_12_ts(ba3,timelen)
ba3_sm = ts_12_sm # historical
##smooth_12_ts(ba3nn,timelen)
##ba3nn_sm = ts_12_sm
lines!(ax, B,ba3_sm[:], 
    linewidth = 0.75,
    color = "lightcyan"
    #label = "GFDL: RONI"
    )

smooth_12_ts(ba4,timelen)
ba4_sm = ts_12_sm # historical
lines!(ax, B,ba4_sm[:], 
    linewidth = 0.75,
    color = "paleturquoise"
    #label = "E3SM: RONI"
    )
#
smooth_12_ts(ba5b,timelen2)
ba5b_sm = ts_12_sm # ssp585
lines!(ax, C,ba5b_sm[:], 
    linewidth = 0.75,
    color = "mistyrose"
    #label = "CESM2: RONI"
    )
#
smooth_12_ts(ba6,timelen)
ba6_sm = ts_12_sm # historical
smooth_12_ts(ba6b,timelen2)
ba6b_sm = ts_12_sm # ssp585
lines!(ax, B,ba6_sm[:], 
    linewidth = 0.75,
    color = "thistle2"
    #label = "HadGEM3: RONI"
    )
lines!(ax, C,ba6b_sm[:], 
    linewidth = 0.75,
    color = "thistle"
    #label = "HadGEM3: RONI"
    )
#
smooth_12_ts(ba9b,timelen2)
ba9b_sm = ts_12_sm # ssp585
lines!(ax, C,ba9b_sm[:], 
    linewidth = 0.75,
    color = "lavender"
    #label = "ACCESS: RONI"
    )
#
#
## ssp585
#tenCent  = [ba1b_sm';ba2b_sm';ba5b_sm';ba9b_sm';ba6b_sm']
# plot individual model for ssp585 if desired.
tenCent  = [ba1b_sm']
mnSSP = ba6b_sm
#mnSSP = mean(tenCent, dims = 1)
## historical
lightOut = [ba1_sm';ba2_sm';ba3_sm';ba4_sm';ba6_sm']
##lightOut = [ba2_sm';ba3_sm';ba4_sm';ba6_sm']
mnH   = mean(lightOut, dims = 1)
#
##smooth_12_ts(ba6,timelen)
##ba6_sm = ts_12_sm # historical
lines!(ax, B,mnH[:], 
    linewidth = 2.0,
    color = "red",
    label = "mn CMIP6"
    )
lines!(ax, C,mnSSP[:], 
    linewidth = 2.0,
    color = "red"
    #label = "mn SSP"
    )
lines!(ax, A,blah3[:], 
    linewidth = 2.,
    color = "black",
    label = "Observed"
    )
#
#
#
##lines!(ax, B,ba3nn_sm[:], 
##    linewidth = 1.5,
##    label = "GFDL: ONI"
##    )
###lines!(ax, nino3p4_anom[:], 
###    linewidth = 1.5,
###    label = "nino 3.4"
###    )
###lines!(ax, ts_oni[:], 
###    linewidth = 1.5,
###    label = "oni"
###    )
###smooth_12_ts(roni_a,2040)
###blah3 = ts_12_sm
####lines!(ax, roni_a[:], 
###lines!(ax, A,blah3[:], 
###    linewidth = 1.5,
###    color = "black",
###    label = "Observed: RONI"
###    )
####axislegend("legend"; position=:rb)
#
axislegend( position=:lt)
#
save("plotENSOcmipHadGEM3.png",fig)








