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

# to run this script on luft:
#/Users/C823281551/.juliaup/bin/julia plotENSOindices.jl

file1 = "/Users/C823281551/data/obs/nina34.noaa.csv"
file2 = "/Users/C823281551/data/obs/oni.noaa.csv"
#file3 = "/Users/C823281551/Downloads/cmip5mean_nino3.4.txt"
# using Excel, I opened the cmip5mean_nino3.4.txt file and then saved it as a csv file.
file3 = "/Users/C823281551/Downloads/cmip5mean_nino3.4test.csv"
file4 = "/Users/C823281551/Downloads/cmip5mean_tropicalmean.csv"

df1 = CSV.read(file1, header = 4, delim="  ", footerskip = 4, DataFrame)
df2 = CSV.read(file2, header = 2, delim="  ", footerskip = 9, DataFrame)
#df3 = CSV.read(file3, header = 0, delim="     ", footerskip = 0, DataFrame)
df3 = CSV.read(file3, header = 0, footerskip = 0, DataFrame) 
df4 = CSV.read(file4, header = 0, footerskip = 0, DataFrame) # seems like you don't need to specify
# 'delim' if the delimiters are comma's.

nms = ["year", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

dfa = DataFrame(df1, nms)
dfb = DataFrame(df2, nms)
dfc = DataFrame(df3, nms)
dfd = DataFrame(df4, nms)

#equals_yr(year::String7) = year == " 2005"

##whichy = [1965, 1970, 1975, 1980, 1985]
#t1 = collect(df1[1, 2:13])
#t2 = collect(df2[1, 2:13])
#t3 = collect(df3[1, 2:13])
#t4 = collect(df4[1, 2:13])

println("~~~~~~~~~")
println(collect(df1[1, 1:13]))
println("~~~~~~~~~")
println(collect(df2[1, 1:13]))
println("~~~~~~~~~")
println(collect(df3[90, 1:13]))
println("~~~~~~~~~")
println(collect(df4[90, 1:13]))
println("~~~~~~~~~")

istart = 2
iend   = 70
for i in istart:iend
    if i < istart + 1 
        global a1 = collect(df1[istart-1, 2:13])
        global a2 = collect(df2[istart-1, 2:13])
        global a3 = collect(df3[istart-1, 2:13])
        global a4 = collect(df4[istart-1, 2:13])
    end
    b1 = collect(df1[i, 2:13]) # grab a row of the DataFrame and convert to vector
    b2 = collect(df2[i, 2:13]) # grab a row of the DataFrame and convert to vector
    b3 = collect(df3[i, 2:13]) # grab a row of the DataFrame and convert to vector
    b4 = collect(df4[i, 2:13]) # grab a row of the DataFrame and convert to vector
    global c1 = [a1; b1] # concatinate two vectors
    global c2 = [a2; b2] # concatinate two vectors
    global c3 = [a3; b3] # concatinate two vectors
    global c4 = [a4; b4] # concatinate two vectors
    a1 = c1
    a2 = c2
    a3 = c3
    a4 = c4
end

# broadcast (using the '.') the parse function to apply to vector.
enso34  = parse.(Float64,c1)
ensooni = parse.(Float64,c2)

ens34c5 = c3 # parse.(Float64,c3) # nino3p4 index from cmip5 models
tmean   = c4 # parse.(Float64,c4)

roni_raw = ens34c5-tmean

print(enso34)
print("*************************")
len = size(enso34)

fig = Figure(;
    size = (700,400),
    )
ax = Axis(fig[1,1];
    xlabel="monthly mean values",
    ylabel="anomaly",
    #xticks=([0:80:800],["1", "2", "3","4","5","6","7","8","9","10"]),
    title="ENSO index comparison"
    )
lines!(ax, enso34[:], 
    linewidth = 2.0,
    label = "Nino 3.4"
    )
lines!(ax, ensooni[:], 
    linewidth = 2.0,
    label = "Oceanic Nino Index"
    )
lines!(ax, roni_raw[:], 
    linewidth = 1.5,
    label = "garbage test"
    )
axislegend("legend"; position=:lt)
#
save("plotENSOinds.png",fig)

