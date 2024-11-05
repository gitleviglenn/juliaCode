#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CMIProni.jl
#
# plot a time series of roni for a particular cmip6 model
#
# levi silvers                                                     nov 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using DataFrames
using CSV
using CairoMakie
using Statistics 
using NCDatasets

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

path="/Users/C823281551/"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Beginning of work with observations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
roni   = zeros(jend)

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

roni = sig_scale.*(ts_oni - tmn_sm)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# beginning of work with GCM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Choose which GCM to use: 

# CNRM
file1  = path*"data/tos_CNRM_hist/tos_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_18500116-20141216.nc"
file1c = path*"data/CNRM-CM6_ssp585_20150116-21001216/tos_Omon_CNRM-CM6-1-HR_ssp585_r1i1p1f2_gn_20150116-21001216remapbil.nc"
file1b = path*"data/CNRM-ESM2_ssp585/tos_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_20150116-21001216.nc"
# MPI-ESM
file1  = path*"data/tos_MPI_hist/tos_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_18500116-20141216_regridded.nc"
#file1b = path*"data/MPI-ESM1_ssp585_20150116-21001216/tos_Omon_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_20150116-21001216_latlon.nc"
# GFDL
#file1b  = path*"data/tos_GFDL_hist/tos_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_18500116-20141216.nc"
# E3SM
#file1  = path*"data/E3SM-1-1-ECA_historical_18500116-20141216/tos_Omon_E3SM-1-1-ECA_historical_r1i1p1f1_gr_18500116-20141216.nc"
#file1b  = path*"data/"
## CESM2
#file1b = path*"data/CESM2_ssp585_20150115-21001215/tos_Omon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215.nc" 
## HadGEM3
#file1  = path*"data/HadGEM3_historical/tos_Omon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_18500116-20141216_regridded.nc"
#file1b = path*"data/HadGEM3_ssp585/tos_Omon_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_20150116-21001216_regridded.nc"
## ACCESS
file1b = path*"data/ACCESS-CM2_ssp585/tos_Omon_ACCESS-CM2_ssp585_r1i1p1f1_gn_20150116-21001216_regridded.nc"

#
timelen = 1980
inpFile = file1
prepare_cmip_ts(inpFile,timelen)
ba1 = ts_rmn
ba1nn = ts_rmn2
#
#MPI scenario timeseries...
timelen2=1032
inpFile = file1c
prepare_cmip_ts(inpFile,timelen2)
ba1c = ts_rmn
inpFile = file1b
prepare_cmip_ts(inpFile,timelen2)
ba1b = ts_rmn

# define thresholds: 
nino_thresh = zeros(timelen+timelen2) .+ 0.75


# make the plot
# define x-axis for time-series plots
A = collect(1854:1/12:2023.92)
B = collect(1850.0833333333333:1/12:2015)
C = collect(2015.083333:1/12:2101)
D = collect(1850.0833333333333:1/12:2101)

fig = Figure(;
    size = (700,400),
    )
ax = Axis(fig[1,1];
    xlabel="monthly mean, smoothed",
    ylabel="ENSO RONI anomalies",
    #xticks=([1850,1860,1870,1880,1890,1900,1910,1920,1930,1940,1950,1960,1970,1980,1990,2000,2010,2020]),
    xticks=([1850,1870,1890,1910,1930,1950,1970,1990,2010,2030,2050,2070,2090]),
    title="ENSO: historical period and ssp585"
    )
smooth_12_ts(roni,2040)
roni_sm = ts_12_sm

# plot the threshold lines
lines!(ax, D,nino_thresh[:],
    color = "red"
    )
lines!(ax, D,-nino_thresh[:],
    color = "red"
    )

# historical
smooth_12_ts(ba1,timelen)
gcm_sm1 = ts_12_sm
lines!(ax, B,gcm_sm1[:], 
    linewidth = 1.25,
    #color = "paleturquoise1"
    color = "black"
    #label = "CNRM: RONI"
    )
limits!(1850, 2100, -4, 4)
# ssp585
smooth_12_ts(ba1b,timelen2)
gcm_sm2 = ts_12_sm
lines!(ax, C,gcm_sm2[:], 
    linewidth = 1.25,
    color = "black"
    #color = "paleturquoise3"
    #label = "CNRM: RONI"
    )

# observed RONI
lines!(ax, A,roni_sm[:], 
    linewidth = 2.,
    color = "brown",
    label = "Observed"
    )

#axislegend( position=:lt)
#
save("plotENSOcmip1modB.png",fig)




