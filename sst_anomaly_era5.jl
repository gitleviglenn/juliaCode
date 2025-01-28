#-----------------------------------------------------------------------------------------------
# sst_anomaly_era5.jl
#
# plot sst from era5
#
# levi silvers                                                                  jan 2025
#-----------------------------------------------------------------------------------------------

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics

path="/Users/C823281551/data/ERA5/"

filein  = path*"era5_sst_1990th2023_360x180.nc"
tag = "ERA5"
data   = NCDataset(filein)

lat = data["lat"]
lon = data["lon"]
tme = data["valid_time"]
 
sst_var = data["sst"]

# create indices to select nino and nina years
#ninoyears = [18 54 90 150 174 234 306 402]
#ninayears = [102 114 210 246 318 366 378 390]
# SH
ninoyears = [23 35 59 96 155 239 311 347]
ninayears = [107 119 215 251 263 335 371 383]
function create_indices(years)
  ensoInd = Matrix{Int64}(undef, 8, 6)
  for i in 1:8
    si = years[i]
    ensoInd[i,:] = si:si+5
  end
  #o1=years[1]:years[1]+5
  #ensoInd = [o1 o2 o3 o4 o5 o6 o7 o8]
  return ensoInd
end 

ninoInd = create_indices(ninoyears)
ninaInd = create_indices(ninayears)

dims = size(sst_var)
numfields = 48
sst_nino          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
sst_nina          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)

lat1 = 50
lat2 = 130

for i in 1:48
  sst_nino[:,lat1:lat2,i] = sst_var[:,lat1:lat2,ninoInd[i]]
  sst_nina[:,lat1:lat2,i] = sst_var[:,lat1:lat2,ninaInd[i]]
end


sst_mn = mean(sst_var, dims=3)
#sst_mn = mean.(skipmissing(sst_var, dims=3))

nino_mn = mean(sst_nino, dims=3)
nina_mn = mean(sst_nina, dims=3)

sst_mn1 = mean(sst_mn, dims=1)
sst_mn2 = mean(sst_mn1, dims=2)

# according to AI, the mean tropical sst is 20C
sst_ai = 293.

# calculate mean between +/-20 degrees: 
sst_tr_mn = mean(skipmissing(sst_var[:,70:110,:]))
nino_composite = nino_mn - nina_mn
nino_anom      = nino_mn[:,:].-sst_mn[:,:,1]
data_2_plot=nina_mn[:,:].-sst_tr_mn
data_3_plot=nino_composite
data_2_plot_tot=sst_mn[:,:,1].-sst_tr_mn

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
#function fig_tot_plot(inpv,inpv2,inpdiff,d1,d2,tit)
function fig_tot_plot(inpv,inpv2,inpv3,d1,d2,tit)
# total sst
# el nino or la nina sst
# nino minus nina or nino minus total
#function fig_tot_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(10,15,10,35),
        backgroundcolor=:white,
        size=(900,1200),
        )
    #ax = Axis(f2[1,1]; #--> default plot is rectangular equidistant 
    #inpdiff = inpv - inpv2
    ax = GeoAxis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        limits=(-180,180,-40,40),
        title="Tropical Mean anom",
        )
        bb = contourf!(ax, d1, d2, inpv, 
             #levels = range(0, 50, length = 25), # tos
             levels = range(-6, 6, length = 21), # rh
             colormap = :vik,
             #colormap = :Blues_8,
             #colormap = :broc,
             #colormap = :bam,
             #colormap = :batlow,
             #colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[1,2], bb)
    ax2 = GeoAxis(f2[2,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        limits=(-180,180,-40,40),
        title="Nino anom",
        )
        bb2 = contourf!(ax2, d1, d2, inpv2, 
             #levels = range(0, 50, length = 25), # tos
             levels = range(-6, 6, length = 21), # rh
             colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        lines!(ax2, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[2,2], bb2)
    ax3 = GeoAxis(f2[3,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        limits=(-180,180,-40,40),
        title="Nino minus Nina",
        )
        bb3 = contourf!(ax3, d1, d2, inpv3, 
             #levels = range(0, 50, length = 25), # tos
             levels = range(-2, 2, length = 21), # rh
             colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        lines!(ax2, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[3,2], bb3)
    return f2
end
tit="ERA5 SST"
fig = fig_anom_plot(nino_composite[:,:,1],lon,lat,tit)
fig1name=tag*"_sst_nino_comp_SH.png"
save(fig1name, fig)
##inpAnom = data_2_plot_tot - data_2_plot
##fig = fig_tot_plot(data_2_plot_tot,data_2_plot,data_3_plot,lon,lat,tit)
fig = fig_tot_plot(data_2_plot_tot,data_2_plot,data_3_plot[:,:,1],lon,lat,tit)
#
#fig = fig_tot_plot(data_2_plot_tot,lon,lat,tit)
fig2name=tag*"_sst_anom_SH.png"
save(fig2name, fig)
