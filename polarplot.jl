#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# polarplot.jl
# 
# plots a 4 panel figure of output from Scythe using a polar plotting projection
# levi silvers                                     august 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
using CairoMakie
using CSV
using DataFrames

# to see colortable options see: 
#https://docs.makie.org/v0.21/explanations/colors

# if running within julia REPL:
# include("/Users/silvers/code/juliaCode/functest.jl")
# if running within a terminal: 
# on luft: 
# /home/C823281551/.julia/juliaup/julia-1.10.4+0.x64.linux.gnu/bin/julia polarplot.jl

path2data = "/bell-scratch/lsilvers/scythe_test/"
exp       = "Twoway_SWslab_wave2/"

fileing1 = path2data*exp*"gridded_out_0.0.csv"
fileing2 = path2data*exp*"gridded_out_16080.0.csv"
fileing3 = path2data*exp*"gridded_out_36000.0.csv"
fileing4 = path2data*exp*"gridded_out_86400.0.csv"

dfg = CSV.read(fileing1, DataFrame)
describe(dfg)
v1 = dfg[:,"w"];
v1a = reshape(v1, (100,601));

dfg = CSV.read(fileing2, DataFrame)
v1 = dfg[:,"w"];
v1b = reshape(v1, (100,601));

dfg = CSV.read(fileing3, DataFrame)
v1 = dfg[:,"w"];
v1c = reshape(v1, (100,601));

dfg = CSV.read(fileing4, DataFrame)
v1 = dfg[:,"w"];
v1d = reshape(v1, (100,601));

size(v1a)

#f = Figure(size = (800, 500))
#ax = PolarAxis(f[1, 1], title = "Surface")
#rs = range(0, 100, 100)
#phis = range(0, 2pi, 601)
##p = surface!(ax, 0..2pi, 0..10, zeros(size(cs)), color = cs, colormap = :coolwarm)
#p = surface!(ax, 0..2pi, 0..100, v1a, colormap = :coolwarm)
##f,p,pltobj = surface!(ax, 0..2pi, 0..100, v1a, colormap = :coolwarm)
#ax.gridz[] = 100
#tightlimits!(ax) # surface plots include padding by default
#Colorbar(f[2, 1], p, vertical = false, flipaxis = false)

#vals=[1, 5, 10, 15, 20, 25, 30, 35, 40]
#vals = range(0,50,41)
vals = range(0,10,41)

f = Figure(size = (1800, 1800))
ax = PolarAxis(f[1, 1], title = "time 1")
#rs = range(0, 100, 100)
rs = range(0, 99, 100)[1:100]
phis = range(0, 2pi, 601)[1:601]
#p = surface!(ax, 0..2pi, 0..100, v1a, colormap = :coolwarm)
p = contourf!(ax, 0..2pi, 0..100, v1a, levels = vals, colormap = :Reds_9)
#p = heatmap!(ax, 0..2pi, 0..100, v1a, colormap = :coolwarm)
#p = voronoiplot!(ax, phis, rs, v1a, show_generators = false, strokewidth = 0, colormap = :coolwarm)
#ax.gridz[] = 100
tightlimits!(ax) # surface plots include padding by default

ax = PolarAxis(f[1, 2], title = "time 2")
p = contourf!(ax, 0..2pi, 0..100, v1b, levels = vals, colormap = :Reds_9)
tightlimits!(ax) # surface plots include padding by default

ax = PolarAxis(f[2, 1], title = "time 3")
p = contourf!(ax, 0..2pi, 0..100, v1c, levels = vals, colormap = :Reds_9)
tightlimits!(ax) # surface plots include padding by default

ax = PolarAxis(f[2, 2], title = "time 4")
p = contourf!(ax, 0..2pi, 0..100, v1d, levels = vals, colormap = :Reds_9)
tightlimits!(ax) # surface plots include padding by default

Colorbar(f[3, 1], p, vertical = false, flipaxis = false)

save("polarplot_w.png", f)

#f
