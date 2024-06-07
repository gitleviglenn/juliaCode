using CairoMakie
using CSV
using DataFrames

# to see colortable options see: 
#https://docs.makie.org/v0.21/explanations/colors

# include("/Users/silvers/code/juliaCode/functest.jl")

fileing1 = "/Users/silvers/code/Scythe.jl/Oneway_SWslab_wave2/gridded_out_120.0.csv"
fileing2 = "/Users/silvers/code/Scythe.jl/Oneway_SWslab_wave2/gridded_out_10200.0.csv"
fileing3 = "/Users/silvers/code/Scythe.jl/Oneway_SWslab_wave2/gridded_out_30000.0.csv"
fileing4 = "/Users/silvers/code/Scythe.jl/Oneway_SWslab_wave2/gridded_out_47160.0.csv"

dfg = CSV.read(fileing1, DataFrame)
describe(dfg)
v1 = dfg[:,"v"];
v1a = reshape(v1, (100,601));

dfg = CSV.read(fileing2, DataFrame)
v1 = dfg[:,"v"];
v1b = reshape(v1, (100,601));

dfg = CSV.read(fileing3, DataFrame)
v1 = dfg[:,"v"];
v1c = reshape(v1, (100,601));

dfg = CSV.read(fileing4, DataFrame)
v1 = dfg[:,"v"];
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
vals = range(0,50,41)

f = Figure(size = (1800, 1800))
ax = PolarAxis(f[1, 1], title = "v1")
#rs = range(0, 100, 100)
rs = range(0, 99, 100)[1:100]
phis = range(0, 2pi, 601)[1:601]
#p = surface!(ax, 0..2pi, 0..100, v1a, colormap = :coolwarm)
p = contourf!(ax, 0..2pi, 0..100, v1a, levels = vals, colormap = :Reds_9)
#p = heatmap!(ax, 0..2pi, 0..100, v1a, colormap = :coolwarm)
#p = voronoiplot!(ax, phis, rs, v1a, show_generators = false, strokewidth = 0, colormap = :coolwarm)
#ax.gridz[] = 100
tightlimits!(ax) # surface plots include padding by default

ax = PolarAxis(f[1, 2], title = "v2")
p = contourf!(ax, 0..2pi, 0..100, v1b, levels = vals, colormap = :Reds_9)
tightlimits!(ax) # surface plots include padding by default

ax = PolarAxis(f[2, 1], title = "v3")
p = contourf!(ax, 0..2pi, 0..100, v1b, levels = vals, colormap = :Reds_9)
tightlimits!(ax) # surface plots include padding by default

ax = PolarAxis(f[2, 2], title = "v4")
p = contourf!(ax, 0..2pi, 0..100, v1b, levels = vals, colormap = :Reds_9)
tightlimits!(ax) # surface plots include padding by default

Colorbar(f[3, 1], p, vertical = false, flipaxis = false)

save("polarplotv4.png", f)

#f
