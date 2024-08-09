#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# grabENSO.jl
#
# - plot an ENSO index verses time
# - initially the nino3.4 index is used 
#
# levi silvers                                              august 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using DataFrames
using CSV
using CairoMakie

fileENSO = "/Users/C823281551/data/obs/nina34.noaa.csv"
#fileENSO = "/Users/C823281551/data/obs/oni.noaa.csv"
dfe = CSV.read(fileENSO, header = 4, delim="  ", footerskip = 4, DataFrame)

nms = ["year", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

dff = DataFrame(dfe, nms)


equals_yr(year::String7) = year == " 2005"

#whichy = [1965, 1970, 1975, 1980, 1985]

istart = 2
iend   = 73
for i in istart:iend
    if i < istart + 1 
        global a = collect(dfe[istart-1, 2:13])
    end
    b = collect(dfe[i, 2:13]) # grab a row of the DataFrame and convert to vector
    global c = [a; b] # concatinate two vectors
    #print(" ",i," ")
    a = c
end

#print(sizeof(c))

# broadcast (using the '.') the parse function to apply to vector.
enso34 = parse.(Float64,c)

print(enso34)
print("*************************")
len = size(enso34)
nino = zeros(len[1])
nina = zeros(len[1])

# define thresholds for the positive and negative phases of ENSO:
up = 1
dn = -1

for i in istart:len[1]
    if enso34[i] > up 
        nino[i] = 1.
    elseif enso34[i] < dn
        nina[i] =-1.
    else
        nino[i] = 0.
        nina[i] = 0.
    end
end

print(nina)

fig = Figure()
ax = Axis(fig[1,1];
    xlabel="monthly mean values",
    ylabel="anomaly",
    title="ENSO index"
    )
lines!(ax, enso34[:], linewidth = 1.5)
lines!(ax, nino[:], linewidth = 3.0)
lines!(ax, nina[:], linewidth = 3.0)
#
save("plotENSOeasy.png",fig)

