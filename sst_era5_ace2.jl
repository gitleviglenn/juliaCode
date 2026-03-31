#-----------------------------------------------------------------------------------------------
# sst_era5_ace2.jl
#
# plot sst from era5 for the years 2005, 2013, and 2024
#
# computes and plots the relative SST (rsst) values, as in Vecchi and Soden's work.
# this means that at each time step (monthly) the tropical mean value, computed between +/- 20 
# of sst is computed and subtracted from each grid point. then, for a particular year, the 
# relative sst field is averaged over time to create 1 map of relative sst.   the anomalous
# sst fields are compputed relative to the relative SST patterns of 2013.  clear as mud? 
# 
# this script also plots the rsst on a regional domain (North Atlantic)
#
# to see a catolog of colormaps: 
# https://juliagraphics.github.io/ColorSchemes.jl/dev/catalogue/
#
# levi silvers                                                             oct 2025
#-----------------------------------------------------------------------------------------------

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics

path="/Users/C823281551/data/ERA5/"

filein   = path*"era5_sst_june2nov_2005_360x180.nc"
filein2  = path*"era5_sst_june2nov_2013_360x180.nc"
filein3  = path*"era5_sst_june2nov_2024_360x180.nc"
tag = "ERA5"
data    = NCDataset(filein)
data2   = NCDataset(filein2)
data3   = NCDataset(filein3)

#lat = data["latitude"]
#lon = data["longitude"]
lat = data["lat"]
lon = data["lon"]
tme = data["valid_time"]
 
sst_var  = data["sst"]
sst_var2 = data2["sst"]
sst_var3 = data3["sst"]

dims = size(sst_var)

lat1 = 50
lat2 = 130

latS = 70
latN = 110

endt = dims[3]

sst_trm            = Array{Union{Missing, Float64}, 1}(undef, endt)
sst_trm2           = Array{Union{Missing, Float64}, 1}(undef, endt)
sst_trm3           = Array{Union{Missing, Float64}, 1}(undef, endt)
sst1               = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], endt)
sst2               = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], endt)
sst3               = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], endt)

for i in 1:endt
  sst_trm[i]  = mean(skipmissing(sst_var[:,latS:latN,i]))
  sst_trm2[i] = mean(skipmissing(sst_var2[:,latS:latN,i]))
  sst_trm3[i] = mean(skipmissing(sst_var3[:,latS:latN,i]))
end

for i in 1:endt
  sst1[:,:,i] = sst_var[:,:,i] .- sst_trm[i]
  sst2[:,:,i] = sst_var2[:,:,i] .- sst_trm2[i]
  sst3[:,:,i] = sst_var3[:,:,i] .- sst_trm3[i]
end

# relative SST values
rsst1_mn  = mean(sst1, dims=3)
rsst2_mn  = mean(sst2, dims=3)
rsst3_mn  = mean(sst3, dims=3)

# SST values (not relative sst)
sst_mn  = mean(sst_var, dims=3)
sst2_mn = mean(sst_var2, dims=3)
sst3_mn = mean(sst_var3, dims=3)

data_2_plot  = sst_mn[:,:].-273.15
data_2_plot2 = sst2_mn[:,:].-273.15
data_2_plot3 = sst3_mn[:,:].-273.15


function fig_anom_reg_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:white,
        size=(600,300),
        )
    #ax = Axis(f2[1,1]; #--> default plot is rectangular equidistant 
    ax = GeoAxis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        limits=(-120,0,0,40),
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             #levels = range(0, 50, length = 25), # tos
             levels = range(-3, 3, length = 21), # rh
             #colormap = :Blues_8,
             #colormap = :broc,
             #colormap = :bam,
             #colormap = :batlow,
             colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[1,2], bb)
    return f2
end
function fig_anom_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:white,
        size=(600,300),
        )
    #ax = Axis(f2[1,1]; #--> default plot is rectangular equidistant 
    ax = GeoAxis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        limits=(-180,180,-40,40),
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             #levels = range(0, 50, length = 25), # tos
             levels = range(-2, 2, length = 21), # rh
             #colormap = :Blues_8,
             #colormap = :broc,
             #colormap = :bam,
             #colormap = :batlow,
             colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[1,2], bb)
    return f2
end
function fig_tot_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:white,
        size=(600,300),
        )
    #ax = Axis(f2[1,1]; #--> default plot is rectangular equidistant 
    ax = GeoAxis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        limits=(-180,180,-40,40),
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             #levels = range(0, 50, length = 25), # tos
             levels = range(10, 30, length = 40), # rh
             #colormap = :Blues_8,
             #colormap = :broc,
             #colormap = :bam,
             colormap = :batlow,
             #colormap = :lajolla,
             #colormap = :matter,
             #colormap = :batlowS,
             #colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[1,2], bb)
    return f2
end
tit="2005 ERA5 SST"
fig = fig_anom_plot(rsst1_mn[:,:],lon,lat,tit)
#fig = fig_tot_plot(data_2_plot[:,:],lon,lat,tit)
fig1name=tag*"_rsst_2005.png"
save(fig1name, fig)
tit="2013 ERA5 SST"
fig = fig_anom_plot(rsst2_mn[:,:],lon,lat,tit)
#fig = fig_tot_plot(data_2_plot2[:,:],lon,lat,tit)
fig1name=tag*"_rsst_2013.png"
save(fig1name, fig)
tit="2024 ERA5 SST"
fig = fig_anom_plot(rsst3_mn[:,:],lon,lat,tit)
#fig = fig_tot_plot(data_2_plot3[:,:],lon,lat,tit)
fig1name=tag*"_rsst_2024.png"
save(fig1name, fig)

tit="2005 ERA5 SST"
rsst1_m_rsst2 = rsst1_mn .- rsst2_mn
fig = fig_anom_reg_plot(rsst1_m_rsst2[:,:],lon,lat,tit)
fig1name=tag*"_rsst_reganom_2005.png"
save(fig1name, fig)

tit="2013 ERA5 SST"
fig = fig_anom_reg_plot(rsst2_mn[:,:],lon,lat,tit)
fig1name=tag*"_rsst_reg_2013.png"
save(fig1name, fig)

rsst3_m_rsst2 = rsst3_mn .- rsst2_mn
tit="2024 ERA5 SST"
fig = fig_anom_reg_plot(rsst3_m_rsst2[:,:],lon,lat,tit)
#fig = fig_anom_reg_plot(rsst3_mn[:,:],lon,lat,tit)
fig1name=tag*"_rsst_reganom_2024.png"
save(fig1name, fig)

data_anom  = data_2_plot .- data_2_plot2
data_anom2 = data_2_plot3 .- data_2_plot2
tit="ERA5 SST anom 2005 - 2013"
fig = fig_anom_plot(data_anom[:,:],lon,lat,tit)
fig1name=tag*"_sst_2005anom.png"
save(fig1name, fig)
tit="ERA5 SST anom 2024 - 2013"
fig = fig_anom_plot(data_anom2[:,:],lon,lat,tit)
fig1name=tag*"_sst_2024anom.png"
save(fig1name, fig)

tit="ERA5 SST 2013"
fig = fig_tot_plot(data_2_plot2[:,:],lon,lat,tit)
fig1name=tag*"_sst_2013.png"
save(fig1name, fig)

