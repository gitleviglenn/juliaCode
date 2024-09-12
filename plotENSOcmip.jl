#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotENSOcmip.jl
#
# - plot ENSO indices for individual cmip models
# - nino 3.4 region is selected
# - anomalies are computed relative to a mean over entire time series.
# - seasonal cycle is removed
# - a three point running mean is applied
#
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

path="/Users/C823281551/"

file1 = path*"data/tos_GFDL_hist/tos_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_18500116-20141216.nc"

ds1 = NCDataset(file1)
ds1.attrib

# to print meta data for a particular variable: 
ds1["tos"]
sst1 = ds1["tos"]
nclat = ds1["lat"]
nclon = ds1["lon"]
nctime = ds1["time"]

# calculate and plot the Nino 3.4 index
# according to psl.noaa.gov, the index is the area average from 5S-5N and 170-120W
# these indices pick out a slightly larger box, by 0.5 on each side:
println(nclat[15:26])
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println(nclon[10:61])
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~")

nino34_full = sst1[10:61, 15:26, :];
typeof(nino34_full)

tlength = 1980
nino34_ts = zeros(tlength)

nino34_ts_mn = mean(skipmissing(nino34_full))

for i = 1:tlength
    nino34_ts[i]=mean(skipmissing(nino34_full[:,:,i]))-nino34_ts_mn
end

# smooth the time series with a running mean
ts_rmn = zeros(tlength)
istart= 2
jend  = tlength
for i in istart:jend-1
  #ts_rmn[i] = (nino34_ts[i+1]+nino34_ts[i]+nino34_ts[i-1])/3
  ts_rmn[i] = nino34_ts[i]
end
ts_rmn[1]=ts_rmn[2]
ts_rmn[jend]=ts_rmn[jend-1]

# remove the seasonal cycle:
# note: these months are not necessarily correct.  it depends on what month the ts starts on...
# 
ss  = zeros(12)
jan = [mean([ts_rmn[i] for i in 1:12:tlength])]
feb = [mean([ts_rmn[i] for i in 2:12:tlength])]
mar = [mean([ts_rmn[i] for i in 3:12:tlength])]
apr = [mean([ts_rmn[i] for i in 4:12:tlength])]
may = [mean([ts_rmn[i] for i in 5:12:tlength])]
jun = [mean([ts_rmn[i] for i in 6:12:tlength])]
jul = [mean([ts_rmn[i] for i in 7:12:tlength])]
aug = [mean([ts_rmn[i] for i in 8:12:tlength])]
sep = [mean([ts_rmn[i] for i in 9:12:tlength])]
oct = [mean([ts_rmn[i] for i in 10:12:tlength])]
nov = [mean([ts_rmn[i] for i in 11:12:tlength])]
dec = [mean([ts_rmn[i] for i in 12:12:tlength])]

ss = [jan feb mar apr may jun jul aug sep oct nov dec]

ts_rmn_nsc = zeros(tlength)
for i in 1:12:tlength
  ts_rmn_nsc[i]    = ts_rmn[i]    - ss[1]
  ts_rmn_nsc[i+1]  = ts_rmn[i+1]  - ss[2]
  ts_rmn_nsc[i+2]  = ts_rmn[i+2]  - ss[3]
  ts_rmn_nsc[i+3]  = ts_rmn[i+3]  - ss[4]
  ts_rmn_nsc[i+4]  = ts_rmn[i+4]  - ss[5]
  ts_rmn_nsc[i+5]  = ts_rmn[i+5]  - ss[6]
  ts_rmn_nsc[i+6]  = ts_rmn[i+6]  - ss[7]
  ts_rmn_nsc[i+7]  = ts_rmn[i+7]  - ss[8]
  ts_rmn_nsc[i+8]  = ts_rmn[i+8]  - ss[9]
  ts_rmn_nsc[i+9]  = ts_rmn[i+9]  - ss[10]
  ts_rmn_nsc[i+10] = ts_rmn[i+10] - ss[11]
  ts_rmn_nsc[i+11] = ts_rmn[i+11] - ss[12]
end

# smooth the time series with a running mean
ts_rmn = zeros(tlength)
istart= 2
jend  = tlength
for i in istart:jend-1
  ts_rmn[i] = (ts_rmn_nsc[i+1]+ts_rmn_nsc[i]+ts_rmn_nsc[i-1])/3
end
ts_rmn[1]=ts_rmn[2]
ts_rmn[jend]=ts_rmn[jend-1]

#fig_plot(nclat,40)
#fig_plot(nino34_ts,tlength,"queen bee")
fig_plot(ts_rmn_nsc,tlength,"Krulle Bol")
#fig_plot(ts_rmn,tlength,"dance floor")





