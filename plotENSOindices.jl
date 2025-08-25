#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotENSOindices.jl
#
# - plot mulitple ENSO indices verses time
#
# levi silvers                                              august 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using DataFrames
using CSV
using CairoMakie
using Statistics
using NCDatasets

# to run this script on luft:
#/Users/C823281551/.juliaup/bin/julia plotENSOindices.jl
# /Users/silvers/.juliaup/bin/julia plotENSOindices.jl

debugSwitch = 0

# include("/Users/silvers/code/juliaCode/plotENSOindices.jl")
#file1 = "/Users/C823281551/data/obs/nina34.noaa.csv"
#file2 = "/Users/C823281551/data/obs/oni.noaa.csv"
#file3 = "/Users/C823281551/Downloads/cmip5mean_nino3.4.txt"
# using Excel, I opened the cmip5mean_nino3.4.txt file and then saved it as a csv file.
# do not save csv files as the UTF-8 version, just the basic .csv
#file3 = "/Users/C823281551/Downloads/cmip5mean_nino3.4test.csv"
#file4 = "/Users/C823281551/Downloads/cmip5mean_tropicalmean.csv"
path="/Users/C823281551/"
#path="/Users/silvers"
file1 = path*"dataB/enso/nino34.noaa.csv"
file2 = path*"dataB/enso/oni.noaa.csv"
file3 = path*"dataB/enso/cmip5_nino3.4.csv"
file4 = path*"dataB/enso/cmip5_tropicalmean.csv"
file5 = path*"dataB/enso/cmip6_nino3.4.csv"
file6 = path*"dataB/enso/cmip6_tropicalmean.csv"
file7 = path*"data/tos_GFDL_hist/tos_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_18500116-20141216.nc"
file8 = path*"data/obs/observed_tropicalmean.csv"


ds1 = Dataset(file7)
ds1.attrib

#sst1 = ds1.tos[:,:,:] # doesn't seem to work for some reason
#
sst1 = ds1["tos"]
nclat = ds1["lat"]
nctime = ds1["time"]
# why doesn't describe(sst1) work?  
println(sst1[1:12])
# sst1[23,10,54]

# typing this in the REPL will output the attributes and meta data:
# ds1.attrib

df1 = CSV.read(file1, header = 4, footerskip = 4, DataFrame)
df2 = CSV.read(file2, header = 2, footerskip = 9, DataFrame)
df8 = CSV.read(file8, header = 0, footerskip = 0, DataFrame)
#df3 = CSV.read(file3, header = 0, delim="     ", footerskip = 0, DataFrame)
df3 = CSV.read(file3, header = 0, footerskip = 0, DataFrame) 
df4 = CSV.read(file4, header = 0, footerskip = 0, DataFrame) # seems like you don't need to specify
df5 = CSV.read(file5, header = 0, footerskip = 0, DataFrame) 
df6 = CSV.read(file6, header = 0, footerskip = 0, DataFrame) # seems like you don't need to specify
# 'delim' if the delimiters are comma's.

nms = ["year", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

dfa  = DataFrame(df1, nms)
dfb  = DataFrame(df2, nms)
dfc  = DataFrame(df3, nms)
dfd  = DataFrame(df4, nms)
dfe  = DataFrame(df5, nms)
dff  = DataFrame(df6, nms)
dftm = DataFrame(df8, nms)

#equals_yr(year::String7) = year == " 2005"

##whichy = [1965, 1970, 1975, 1980, 1985]
#t1 = collect(df1[1, 2:13])
#t2 = collect(df2[1, 2:13])
#t3 = collect(df3[1, 2:13])
#t4 = collect(df4[1, 2:13])

if (debugSwitch < 1)
    println("~~~~~~~~~")
    println(collect(df1[1, 1:13]))
    println("~~~~~~~~~")
    println(collect(df2[1, 1:13]))
    println("~~~observed trop mn data: ~~~~")
    println(collect(df8[1, 1:13]))
    println("~~~~~~~~~")
    println(collect(df3[91, 1:13]))
    println("~~~~~~~~~")
    println(collect(df4[91, 1:13]))
    println("~~~~~~~~~")
    println(collect(df5[102, 1:13]))
    println("~~~~~~~~~")
    println(collect(df6[102, 1:13]))
    println("~~~~cmip5~~~~~")
    println(collect(df3[130, 1:13]))
    println(collect(df3[160, 1:13]))
    println("~~~~cmip5~~~~~")
    println(collect(df4[130, 1:13]))
    println(collect(df4[160, 1:13]))
    println("~~~~cmip6~~~~~")
    println(collect(df5[141, 1:13]))
    println(collect(df5[171, 1:13]))
    println("~~~~cmip6~~~~~")
    println(collect(df6[141, 1:13]))
    println(collect(df6[171, 1:13]))
    println("~~~~~~~~~")
end

istart = 2
iend   = 70
for i in istart:iend
    if i < istart + 1 
        global a1 = collect(df1[istart-1, 2:13])
        global a2 = collect(df2[istart-1, 2:13])
    end
    b1 = collect(df1[i, 2:13]) # grab a row of the DataFrame and convert to vector
    b2 = collect(df2[i, 2:13]) # grab a row of the DataFrame and convert to vector
    global c1 = [a1; b1] # concatinate two vectors
    global c2 = [a2; b2] # concatinate two vectors
    a1 = c1
    a2 = c2
end

# cmip5 data
istart = 92
iend   = 239
for i in istart:iend
    if i < istart + 1 
        global a3 = collect(df3[istart-1, 2:13])
        global a4 = collect(df4[istart-1, 2:13])
    end
    b3 = collect(df3[i, 2:13]) # grab a row of the DataFrame and convert to vector
    b4 = collect(df4[i, 2:13]) # grab a row of the DataFrame and convert to vector
    global c3 = [a3; b3] # concatinate two vectors
    global c4 = [a4; b4] # concatinate two vectors
    a3 = c3
    a4 = c4
end

# cmip6 data
istart = 103
iend   = 250
for i in istart:iend
    if i < istart + 1 
        global a5 = collect(df5[istart-1, 2:13])
        global a6 = collect(df6[istart-1, 2:13])
    end
    b5 = collect(df5[i, 2:13]) # grab a row of the DataFrame and convert to vector
    b6 = collect(df6[i, 2:13]) # grab a row of the DataFrame and convert to vector
    global c5 = [a5; b5] # concatinate two vectors
    global c6 = [a6; b6] # concatinate two vectors
    a5 = c5
    a6 = c6
end

# subtract this mean to get anomalies
basec3=c3[130:160]
basec4=c4[130:160]
basec5=c5[141:171]
basec6=c6[141:171]
mnc3 = mean(basec3)
mnc4 = mean(basec4)
mnc5 = mean(basec5)
mnc6 = mean(basec6)

# anomalies
c3 = c3 .- mnc3
c4 = c4 .- mnc4
c5 = c5 .- mnc5
c6 = c6 .- mnc6

# std (before or after removing seasonal cycle?)
sig_oni = std(c3)
sig_diff = std(c3-c4)
sig_scale = sig_oni/sig_diff

jend = 148*12+12
jend2 = 148*12+12
#jend=iend*12

# create arrays for the running mean time series
rmn_cm5 = zeros(jend)
rmn_cm6 = zeros(jend)

# describe(a3)
# describe(rmn)

## calculate the 3 month running mean of nino34
#istart = 2
#for i in istart:jend-1
#  rmn_cm5[i] = (c3[i+1]+c3[i]+c3[i-1])/3
#  rmn_cm6[i] = (c5[i+1]+c5[i]+c5[i-1])/3
#end
#rmn_cm5[1]=rmn_cm5[2]
#rmn_cm5[jend]=rmn_cm5[jend-1]
#rmn_cm6[1]=rmn_cm6[2]
#rmn_cm6[jend]=rmn_cm6[jend-1]

# get the monthly mean values to remove the seasonal cycle: 
# df3 --> nino3p4 for cmip5 models
# df4 --> tropical mean for cmip5 models
c3nsc  = zeros(jend)
c4nsc  = zeros(jend)
c5nsc  = zeros(jend)
c6nsc  = zeros(jend)
#tr_nsc = zeros(jend)
mmn3 = zeros(12)
mmn3 = [mean(df3[:,i]) for i in 2:13]
mmn = zeros(12)
mmn = [mean(df4[:,i]) for i in 2:13]
mmn_34_cm6 = zeros(12)
mmn_34_cm6 = [mean(df5[:,i]) for i in 2:13]
mmn_cm6 = zeros(12)
mmn_cm6 = [mean(df6[:,i]) for i in 2:13]
# compute seasonal cycle for the observed tropical mean:
mmn_tr = zeros(12)
mmn_tr = [mean(df8[:,i]) for i in 2:13]
#
#endpt = 148*12+12
# remove the seasonal cycle
for i in 1:12:jend
  # cmip5
  c3nsc[i]=c3[i]-mmn3[1]
  c3nsc[i+1]=c3[i+1]-mmn3[2]
  c3nsc[i+2]=c3[i+2]-mmn3[3]
  c3nsc[i+3]=c3[i+3]-mmn3[4]
  c3nsc[i+4]=c3[i+4]-mmn3[5]
  c3nsc[i+5]=c3[i+5]-mmn3[6]
  c3nsc[i+6]=c3[i+6]-mmn3[7]
  c3nsc[i+7]=c3[i+7]-mmn3[8]
  c3nsc[i+8]=c3[i+8]-mmn3[9]
  c3nsc[i+9]=c3[i+9]-mmn3[10]
  c3nsc[i+10]=c3[i+10]-mmn3[11]
  c3nsc[i+11]=c3[i+11]-mmn3[12]
  c4nsc[i]=c4[i]-mmn[1]
  c4nsc[i+1]=c4[i+1]-mmn[2]
  c4nsc[i+2]=c4[i+2]-mmn[3]
  c4nsc[i+3]=c4[i+3]-mmn[4]
  c4nsc[i+4]=c4[i+4]-mmn[5]
  c4nsc[i+5]=c4[i+5]-mmn[6]
  c4nsc[i+6]=c4[i+6]-mmn[7]
  c4nsc[i+7]=c4[i+7]-mmn[8]
  c4nsc[i+8]=c4[i+8]-mmn[9]
  c4nsc[i+9]=c4[i+9]-mmn[10]
  c4nsc[i+10]=c4[i+10]-mmn[11]
  c4nsc[i+11]=c4[i+11]-mmn[12]
end
for i in 1:12:jend2
  # cmip6
  c5nsc[i]=c5[i]-mmn_34_cm6[1]
  c5nsc[i+1]=c5[i+1]-mmn_34_cm6[2]
  c5nsc[i+2]=c5[i+2]-mmn_34_cm6[3]
  c5nsc[i+3]=c5[i+3]-mmn_34_cm6[4]
  c5nsc[i+4]=c5[i+4]-mmn_34_cm6[5]
  c5nsc[i+5]=c5[i+5]-mmn_34_cm6[6]
  c5nsc[i+6]=c5[i+6]-mmn_34_cm6[7]
  c5nsc[i+7]=c5[i+7]-mmn_34_cm6[8]
  c5nsc[i+8]=c5[i+8]-mmn_34_cm6[9]
  c5nsc[i+9]=c5[i+9]-mmn_34_cm6[10]
  c5nsc[i+10]=c5[i+10]-mmn_34_cm6[11]
  c5nsc[i+11]=c5[i+11]-mmn_34_cm6[12]
  c6nsc[i]=c6[i]-mmn_cm6[1]
  c6nsc[i+1]=c6[i+1]-mmn_cm6[2]
  c6nsc[i+2]=c6[i+2]-mmn_cm6[3]
  c6nsc[i+3]=c6[i+3]-mmn_cm6[4]
  c6nsc[i+4]=c6[i+4]-mmn_cm6[5]
  c6nsc[i+5]=c6[i+5]-mmn_cm6[6]
  c6nsc[i+6]=c6[i+6]-mmn_cm6[7]
  c6nsc[i+7]=c6[i+7]-mmn_cm6[8]
  c6nsc[i+8]=c6[i+8]-mmn_cm6[9]
  c6nsc[i+9]=c6[i+9]-mmn_cm6[10]
  c6nsc[i+10]=c6[i+10]-mmn_cm6[11]
  c6nsc[i+11]=c6[i+11]-mmn_cm6[12]
end

sig_oniB = std(c3nsc)
sig_diffB = std(c3nsc-c4nsc)
sig_scaleB = sig_oniB/sig_diffB

sig_oni_cm6   = std(c5nsc)
sig_diff_cm6  = std(c5nsc-c6nsc)
sig_scale_cm6 = sig_oni_cm6/sig_diff_cm6

#for i in istart:jend-1
istart= 2
jend  = 1788
for i in istart:jend-1
  rmn_cm5[i] = (c3nsc[i+1]+c3nsc[i]+c3nsc[i-1])/3
  rmn_cm6[i] = (c5nsc[i+1]+c5nsc[i]+c5nsc[i-1])/3
end
rmn_cm5[1]=rmn_cm5[2]
rmn_cm6[1]=rmn_cm6[2]
rmn_cm5[jend]=rmn_cm5[jend-1]
rmn_cm6[jend]=rmn_cm6[jend-1]

# broadcast (using the '.') the parse function to apply to vector.
#enso34  = parse.(Float64,c1)
#ensooni = parse.(Float64,c2)
enso34  = c1
ensooni = c2

ens34c5 = c3 # parse.(Float64,c3) # nino3p4 index from cmip5 models
# tropical mean (what latitudes?  was land included?)
tmean   = c4 # parse.(Float64,c4)
ens34_cm6 = c5
tmean_cm6 = c6

#roni_raw = ens34c5-tmean
#roni_t   = c3nsc-c4nsc
#roni_t   = rmn-tmean
roni_t   = sig_scaleB.*(rmn_cm5-c4nsc)
roni_cm6   = sig_scale_cm6.*(rmn_cm6-c6nsc)

#print(enso34)
print("*************************")
len = size(enso34)

fig = Figure(;
    size = (700,400),
    )
ax = Axis(fig[1,1];
    xlabel="monthly mean values",
    ylabel="anomaly",
    #xticks=([108,228,348,468,588,708,828,948],["1960","1970","1980","1990","2000","2010","2020","2030"]),
    xticks=([108,228,348,468,588,708,828,948,1068,1188,1308,1428,1548],["1960","1970","1980","1990","2000","2010","2020","2030","2040","2050","2060","2070","2080"]),
    title="ENSO index comparison"
    )
limits!(0, 1700, -4, 7)
lines!(ax, enso34[:], 
    linewidth = 2.0,
    label = "Nino 3.4"
    )
lines!(ax, ensooni[:], 
    linewidth = 2.0,
    label = "Oceanic Nino Index"
    )
lines!(ax, roni_t[:], 
    linewidth = 1.5,
    label = "RONI cmip5 tas"
    )
lines!(ax, roni_cm6[:], 
    linewidth = 1.5,
    label = "RONI cmip6 tas"
    )
axislegend("legend"; position=:rb)
#
save("plotENSOinds_cmip6_1.png",fig)

