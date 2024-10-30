#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# tcRIearthMap.jl
#
# - plot the points where RI first occurs on a world map
# - specifically, read a csv file, use DataFrames to manipulate and plot
# the data.
# - use GeoAxis to plot data on a physical map of the Earth between +/-40
# - particular basins can be highlighted if desired
#
# levi silvers                                              august 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
using DataFrames
using CSV
using CairoMakie
using LaTeXStrings
using GeoMakie

# see various projection options here: 
# https://proj.org/en/9.4/operations/projections/index.html

# run with: 
# /Users/silvers/.juliaup/bin/julia filename.jl

# /Users/silvers/data/hurricane_data/AtlanticTCACE18512023.csv 

# incoming file: 
#filein = "/Users/silvers/data/hurricane_data/FirstRILocationbyStorm.csv"
filein = "/Users/C823281551/data/hurricane_data/FirstRILocationbyStorm.csv"
file2in = "/Users/C823281551/data/hurricane_data/RI_Episodes.csv"


# read a csv file into a dataframe:
#df1 = DataFrame(CSV.File(filein))

dfRI  = CSV.read(filein, DataFrame)
dfRIe = CSV.read(file2in, DataFrame)

# look at the metadata for dfRImap :
describe(dfRI)
# one can also describe  particular columns: 
#describe(dfRI, cols=3:10)

img = rotr90(GeoMakie.earth())

fig = Figure(
    size = (1920, 880),
)
#ga = GeoAxis(
#    fig[1, 1]; # any cell of the figure's layout
#    dest = "+proj=wintri", # the CRS in which you want to plot
#)
ga = GeoAxis(
    fig[1, 1]; # any cell of the figure's layout
    source="+proj=latlong",
    dest = "+proj=eqearth", # the CRS in which you want to plot
    xticks = -180:10:180,
    yticks = -90:10:90,
)
image!(ga, -180..180, -90..90, img; interpolate = false) 
#lines!(ga, GeoMakie.coastlines()) # plot coastlines from Natural Earth as a reference

# You can plot your data the same way you would in Makie
#scatter!(ga, -120:15:120, -60:7.5:60; color = -60:7.5:60, strokecolor = (:black, 0.2))

# looking at these storms (763-1110) corresponds roughly to the events between 2014-2014.
p1=1
#p2=1110
p2=4490

#xRI = dfRI[p1:p2, :LATITUDE]
#yRI = dfRI[p1:p2, :LONGITUDE]
xRI = dfRIe[p1:p2, :USA_LAT]
yRI = dfRIe[p1:p2, :USA_LON]
zRI = dfRIe[p1:p2, :"24-Hr Intensification Rate"]
yearRI = dfRI[:, :YEAR]
#scatter(yRI,xRI)
#scatter!(ga, yRI,xRI, markersize = 7, color = :black)
#scatter!(ga, yRI,xRI,color = dfRIe.var"24-Hr Intensification Rate", markersize = 7, colormap =:batlow) # --> still need to couple the magnitude with the 
#scatter!(ga, yRI,xRI, color=zRI, markersize = 9, colormap =:batlow) # --> still need to couple the magnitude with the 
scatter!(ga, yRI,xRI, color=zRI, markersize = 9, colormap =:heat) # --> still need to couple the magnitude with the 
# color somehow.   dfRIe.var"24-Hr Intensification Rate"

ylims!(ga, -30, 40)
#xlims!(ga, -100, -70)

## select by basin: 
#dfNA=dfRI[in.(dfRI.BASIN, Ref(["NA"])), :]
#xRI = dfNA[:, :LATITUDE]
#yRI = dfNA[:, :LONGITUDE]
#yearRI = dfNA[:, :YEAR]
#scatter!(ga, yRI,xRI, markersize = 7, color = :red)

save("tcRIearthMapE.png", fig)

