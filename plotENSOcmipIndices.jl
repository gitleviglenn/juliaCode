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
# - to use the data from MPI-ESM some regridding was necessary, using both
# - nco and cdo functions: 
# ncks -d time,0,1979 tos_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_18500116-20141216.nc test_1980ts.nc
# ncatted -a coordinates,tos,c,c,"latitude longitude" test_1980ts.nc
# cdo -L remapbil,mygrid -sethalo,-1,-1 test_1980ts.nc test_1980ts_latlon3.nc
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
# functions:  fig_plot(), prepare_cmip_ts(), smooth_ts()
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
        size=(1000,300),
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

function prepare_cmip_ts(inpFile,len)
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
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~lat points: ~~~~~~~~~~~~~~~~")
    println(nclat[15:26])
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
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
    
    nino34_ts_mn = mean(skipmissing(nino34_full))
    tropmn_ts_mn = mean(filter(!isnan, skipmissing(trop_full)))
    
    #println("fury fury fury")
    #println(tropmn_ts_mn) 
    # compute the anomalies
    for i = 1:tlength
        nino34_ts[i]=mean(skipmissing(nino34_full[:,:,i]))-nino34_ts_mn
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

#file1 = path*"data/tos_GFDL_hist/tos_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_18500116-20141216.nc"
#file1 = path*"data/tos_MPI_hist/tos_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_18500116-20141216.nc"
file1 = path*"data/tos_CNRM_hist/tos_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_18500116-20141216.nc"
file2 = path*"data/tos_MPI_hist/tos_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_18500116-20141216_regridded.nc"
file2b = path*"data/MPI-ESM1_ssp585_20150116-21001216/tos_Omon_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_20150116-21001216_latlon.nc"
file3 = path*"data/tos_GFDL_hist/tos_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_18500116-20141216.nc"

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
prepare_cmip_ts(inpFile,timelen)
ba1 = ts_rmn
ba1nn = ts_rmn2
#
inpFile = file2
prepare_cmip_ts(inpFile,timelen)
ba2 = ts_rmn
ba2nn = ts_rmn2

inpFile = file3
prepare_cmip_ts(inpFile,timelen)
ba3 = ts_rmn
ba3nn = ts_rmn2

#MPI scenario timeseries...
timelen2=1032
inpFile = file2b
prepare_cmip_ts(inpFile,timelen2)
ba2b = ts_rmn

timelen = 1980
inpFile = file3
prepare_cmip_ts(inpFile,timelen)
ba3 = ts_rmn

A = collect(1854:1/12:2023.92)
B = collect(1850.0833333333333:1/12:2015)

fig = Figure(;
    size = (700,400),
    )
ax = Axis(fig[1,1];
    xlabel="monthly mean, smoothed",
    ylabel="ENSO index anomalies",
    xticks=([1850,1860,1870,1880,1890,1900,1910,1920,1930,1940,1950,1960,1970,1980,1990,2000,2010,2020]),
    title="ENSO over the historical period"
    )
smooth_12_ts(roni_a,2040)
blah3 = ts_12_sm
lines!(ax, A,blah3[:], 
    linewidth = 2.5,
    color = "black",
    label = "Observed: RONI"
    )
limits!(1850, 2030, -4, 4)

smooth_12_ts(ba1,timelen)
ba1_sm = ts_12_sm
smooth_12_ts(ba1nn,timelen)
ba1nn_sm = ts_12_sm
lines!(ax, B,ba1_sm[:], 
    linewidth = 1.5,
    label = "CNRM: RONI"
    )

smooth_12_ts(ba2,timelen)
ba2_sm = ts_12_sm
smooth_12_ts(ba2nn,timelen)
ba2nn_sm = ts_12_sm
lines!(ax, B,ba2_sm[:], 
    linewidth = 1.5,
    label = "MPI: RONI"
    )
smooth_12_ts(ba3,timelen)
ba3_sm = ts_12_sm
smooth_12_ts(ba3nn,timelen)
ba3nn_sm = ts_12_sm
lines!(ax, B,ba3_sm[:], 
    linewidth = 1.5,
    label = "GFDL: RONI"
    )
#lines!(ax, B,ba3nn_sm[:], 
#    linewidth = 1.5,
#    label = "GFDL: ONI"
#    )
##lines!(ax, nino3p4_anom[:], 
##    linewidth = 1.5,
##    label = "nino 3.4"
##    )
##lines!(ax, ts_oni[:], 
##    linewidth = 1.5,
##    label = "oni"
##    )
##smooth_12_ts(roni_a,2040)
##blah3 = ts_12_sm
###lines!(ax, roni_a[:], 
##lines!(ax, A,blah3[:], 
##    linewidth = 1.5,
##    color = "black",
##    label = "Observed: RONI"
##    )
###axislegend("legend"; position=:rb)
axislegend( position=:lt)
#
save("plotENSOcmipInds.png",fig)








