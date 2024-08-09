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

file1 = "/Users/C823281551/data/obs/nina34.noaa.csv"
file2 = "/Users/C823281551/data/obs/oni.noaa.csv"

df1 = CSV.read(file1, header = 4, delim="  ", footerskip = 4, DataFrame)
df2 = CSV.read(file2, header = 2, delim="  ", footerskip = 9, DataFrame)

nms = ["year", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

dfa = DataFrame(df1, nms)
dfb = DataFrame(df2, nms)

#equals_yr(year::String7) = year == " 2005"

#whichy = [1965, 1970, 1975, 1980, 1985]
t1 = collect(df1[1, 2:13])
t2 = collect(df2[1, 2:13])

println("~~~~~~~~~")
println(collect(df1[1, 1:13]))
println("~~~~~~~~~")
println(collect(df2[1, 1:13]))
println("~~~~~~~~~")

istart = 2
iend   = 70
for i in istart:iend
    if i < istart + 1 
        global a1 = collect(df1[istart-1, 2:13])
        global a2 = collect(df1[istart-1, 2:13])
    end
    b1 = collect(df1[i, 2:13]) # grab a row of the DataFrame and convert to vector
    b2 = collect(df2[i, 2:13]) # grab a row of the DataFrame and convert to vector
    global c1 = [a1; b1] # concatinate two vectors
    global c2 = [a2; b2] # concatinate two vectors
    a1 = c1
    a2 = c2
end

# broadcast (using the '.') the parse function to apply to vector.
enso34  = parse.(Float64,c1)
ensooni = parse.(Float64,c2)

print(enso34)
print("*************************")
len = size(enso34)

fig = Figure()
ax = Axis(fig[1,1];
    xlabel="monthly mean values",
    ylabel="anomaly",
    title="ENSO index comparison"
    )
lines!(ax, enso34[:], linewidth = 2.0)
lines!(ax, ensooni[:], linewidth = 2.0)
#
save("plotENSOinds.png",fig)

