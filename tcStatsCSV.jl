using DataFrames
using CSV
using CairoMakie
using LaTeXStrings
#using Plots
#using Statistics

# run with: 
# /Users/silvers/.juliaup/bin/julia tcStatsCSV.jl

# /Users/silvers/data/hurricane_data/AtlanticTCACE18512023.csv 

# incoming file: 
#filein = "/Users/silvers/data/hurricane_data/AtlanticTCACE18512023.csv"
filein = "/Users/C823281551/data/hurricane_data/AtlanticTCACE18512023.csv"

# read a csv file into a dataframe:
#df1 = DataFrame(CSV.File(filein))

df2 = CSV.read(filein, DataFrame)

# how do I look at the metadata for df?
describe(df2)
# one can also describe  particular columns: 
describe(df2, cols=3:10)

# create a scatter plot:
#Axis(fig[1,1],
fig = Figure()
    ax = Axis(fig[1,1];
        xlabel="days of the year",
        ylabel=latexstring("10^4")*" knots"*latexstring("^2"),
        title="Atlantic ACE",
        )
  #title!(ax1,"in the aeroplane over the sea")
  #xlims!(ax1,0, 365)
  #ylims!(ax1,0, 200)
  #xlabel!(ax1,"up and over")
  #ylabel!(ax1,"bro")
#)

#CairoMakie.scatter(df2[:, "2000"], markersize = 5.0)
#fig, ax, sc = scatter(df2[:, "2005"], linewidth = 5.0)
        lines!(ax, df2[:, "2005"], linewidth = 2.0)
#fig, ax, sc = plot(0:370, [df2[:, "2005"] df2[:, "1915"]])
#xlims!(ax, 0, 370)
        lines!(ax, df2[:, "2006"], linewidth = 2.5)
        lines!(ax, df2[:, "2007"], linewidth = 3.0)
        lines!(ax, df2[:, "2008"], linewidth = 3.5)
        lines!(ax, df2[:, "2009"], linewidth = 4.0)
        lines!(ax, df2[:, "2010"], linewidth = 4.5)
        lines!(ax, df2[:, "2011"], linewidth = 5.0)
        lines!(ax, df2[:, "2012"], linewidth = 5.5)
        lines!(ax, df2[:, "2013"], linewidth = 6.0)
        lines!(ax, df2[:, "2014"], linewidth = 6.5)
        lines!(ax, df2[:, "2015"], linewidth = 7.0)

#x = range(1, 365, length = 365)
#y1 = df2[:, "1917"]
#y2 = df2[:, "1918"]
#y3 = df2[:, "1919"]
#fig, ax, sc = CairoMakie.plot(x, [y1 y2 y3])
#xlims!(ax, 180, 365)

# what is the difference between using CSV.File and CSV.read?  

save("lineplot.png", fig)

# this creates an interesting heatmapsort of plot: 
#CairoMakie.plot([df2[:, "2005"] df2[:, "1915"] df2[:, "1900"]])
# but I am struggling tomake it work for all the years.
