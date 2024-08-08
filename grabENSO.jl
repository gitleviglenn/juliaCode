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

fileENSO = "/Users/C823281551/nina34.anom.data"
dfe = CSV.read(fileENSO, header = 4, delim="  ", footerskip = 4, DataFrame)

nms = ["year", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

dff = DataFrame(dfe, nms)


equals_yr(year::String7) = year == " 2005"

whichy = [1965, 1970, 1975, 1980, 1985]

istart = 63
iend   = 73
for i in istart:73
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

fig = Figure()
ax = Axis(fig[1,1];
    xlabel="monthly mean values",
    ylabel="anomaly",
    title="ENSO index over 12 years"
    )
lines!(ax, enso34[:], linewidth = 3.0)
#
save("plotENSOeasy.png",fig)

