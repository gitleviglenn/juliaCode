using DataFrames
using CSV
using CairoMakie
using Plots

# run with
# /Users/silvers/.juliaup/bin/julia tcRIhist.jl

filein = "/Users/silvers/data/hurricane_data/FirstRILocationbyStorm.csv"

dfRI = CSV.read(filein, DataFrame)  
describe(dfRI) 

fig = Figure()
ax = Axis(fig[1,1];
    xlabel = "year",
    ylabel = "# of storms",
    title  = "Storms undergoing Rapid Intensification",
    )

#histogram!(ax, dfRI.YEAR, label = false, bins=34, color=:gray)
#histogram!(ax, dfRI.YEAR)
#plot(histogram!(ax, dfRI.YEAR, label = false, bins = 34, color=:gray))

# I think 'hist' is from Makie, not Plots.  Not sure why hist works and 
# not histogram.   whatever you schmuck.  
hist!(ax, dfRI.YEAR, label = false, bins=34, color=:gray)


save("tcRIhist.png", fig)

# dfNI=dfRI[in.(dfRI.BASIN, Ref(["NI"])), :]
