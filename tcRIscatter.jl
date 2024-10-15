# plot the points where RI first occurs
using DataFrames
using CSV
using CairoMakie
using LaTeXStrings
using GeoMakie
#using Statistics

# see various projection options here: 
# https://proj.org/en/9.4/operations/projections/index.html

# run with: 
# /Users/silvers/.juliaup/bin/julia tcStatsCSV.jl

# /Users/silvers/data/hurricane_data/AtlanticTCACE18512023.csv 

# incoming file: 
#filein = "/Users/silvers/data/hurricane_data/FirstRILocationbyStorm.csv"
filein = "/Users/C823281551/data/hurricane_data/FirstRILocationbyStorm.csv"

# read a csv file into a dataframe:
#df1 = DataFrame(CSV.File(filein))

dfRI = CSV.read(filein, DataFrame)

# how do I look at the metadata for df?
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
)
image!(ga, -180..180, -90..90, img; interpolate = false) 
#image!(ga, -180:10:180, -90:10:90, img; interpolate = false) 
lines!(ga, GeoMakie.coastlines()) # plot coastlines from Natural Earth as a reference

# You can plot your data the same way you would in Makie
#scatter!(ga, -120:15:120, -60:7.5:60; color = -60:7.5:60, strokecolor = (:black, 0.2))

# looking at these storms (763-1110) corresponds roughly to the events between 2014-2014.
p1=1
p2=1110

xRI = dfRI[p1:p2, :LATITUDE]
yRI = dfRI[p1:p2, :LONGITUDE]
yearRI = dfRI[:, :YEAR]
#scatter(yRI,xRI)
scatter!(ga, yRI,xRI, markersize = 7, color = :black)

ylims!(ga, -40, 40)

# select by basin: 
dfNA=dfRI[in.(dfRI.BASIN, Ref(["NA"])), :]
xRI = dfNA[:, :LATITUDE]
yRI = dfNA[:, :LONGITUDE]
yearRI = dfNA[:, :YEAR]
scatter!(ga, yRI,xRI, markersize = 7, color = :black)


save("tcRIscatterBlackB.png", fig)

