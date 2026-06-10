#-----------------------------------------------------------------------------------------------
# test_mapHR_era5.jl
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
# in the Atlantic basin, the main development region is defined to be from 
# 10-20N and 85-20W.
#
# it looks like the first longitude point is at the prime meridian.
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

#filein   = path*"era5_sst_june2nov_2005_360x180.nc"
#filein2  = path*"era5_sst_june2nov_2013_360x180.nc"
#filein3  = path*"era5_sst_june2nov_2024_360x180.nc"
filein   = path*"era5_sst_june2nov_2005.nc"
filein2  = path*"era5_sst_june2nov_2013.nc"
filein3  = path*"era5_sst_june2nov_2024.nc"
tag = "ERA5"
# lr for low resolution; hr for high resolution
res_tag="_hr"
data    = NCDataset(filein)
data2   = NCDataset(filein2)
data3   = NCDataset(filein3)

lat = data["latitude"]
lon_orig = data["longitude"]
#lat = data["lat"]
#lon = data["lon"]
tme = data["valid_time"]
 
sst_var  = data["sst"]
sst_var2 = data2["sst"]
sst_var3 = data3["sst"]

println("attributes of sst_var are: ",sst_var.attrib)

# Convert longitude from 0-360 to -180-180 convention
lon = Float64.(Array(lon_orig))
lon[lon .> 180] .-= 360
# Find indices to reorder (0-180 goes to end, 180-360 goes to start)
split_idx = findfirst(x -> x > 180, Array(lon_orig))
reorder_idx = vcat(split_idx:length(lon), 1:(split_idx-1))

# Reorder data and longitude
sst_var = sst_var[reorder_idx, :, :]
sst_var2 = sst_var2[reorder_idx, :, :]
sst_var3 = sst_var3[reorder_idx, :, :]
lon = lon[reorder_idx]
lon = sort(lon)

dims = size(sst_var)

println("dims of sst_var are: ",dims)

println("whare the selected latitudes? ",lat[:])

lat1 = 50
lat2 = 130

latS = 70
latN = 110
#latS = findfirst(x -> x == latS_int,lat)
#latN = findfirst(x -> x == latN_int,lat)

println("whare the selected latitudes? ",lat[latS:latN])

println("dims of lat are: ",size(lat))
println("central lat is: ",lat[359:361])
println("tropics? lat is: ",lat[280:440])
latS = 281
latN = 441

endt = dims[3]
println("number of time steps is: ",endt)

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

# mean SST in Celcius
data_2_plot  = sst_mn[:,:].-273.15
data_2_plot2 = sst2_mn[:,:].-273.15
data_2_plot3 = sst3_mn[:,:].-273.15

# Regional average computation
lat_min = 10
lat_max = 20
lon_min = -85
lon_max = -20

# Find indices corresponding to the region
lat_idx = findall(x -> x >= lat_min && x <= lat_max, Array(lat))
lon_idx = findall(x -> x >= lon_min && x <= lon_max, Array(lon))

# Compute regional average for each year
regional_avg_2005 = mean(skipmissing(sst_var[lon_idx, lat_idx, :]))
regional_avg_2013 = mean(skipmissing(sst_var2[lon_idx, lat_idx, :]))
regional_avg_2024 = mean(skipmissing(sst_var3[lon_idx, lat_idx, :]))

# Convert to Celsius
regional_avg_2005_c = regional_avg_2005 - 273.15
regional_avg_2013_c = regional_avg_2013 - 273.15
regional_avg_2024_c = regional_avg_2024 - 273.15
# Compute regional average of relative SST for each year
regional_rsst_avg_2005 = mean(skipmissing(sst1[lon_idx, lat_idx, :]))
regional_rsst_avg_2013 = mean(skipmissing(sst2[lon_idx, lat_idx, :]))
regional_rsst_avg_2024 = mean(skipmissing(sst3[lon_idx, lat_idx, :]))

println("Regional Average SST ($(lat_min)°N to $(lat_max)°N, $(lon_min)°E to $(lon_max)°E):")
println("2005: $(round(regional_avg_2005_c, digits=2)) °C")
println("2013: $(round(regional_avg_2013_c, digits=2)) °C")
println("2024: $(round(regional_avg_2024_c, digits=2)) °C")
println("Regional Average of Relative SST ($(lat_min)°N to $(lat_max)°N, $(lon_min)°E to $(lon_max)°E):")
println("2005: $(round(regional_rsst_avg_2005, digits=2)) K")
println("2013: $(round(regional_rsst_avg_2013, digits=2)) K")
println("2024: $(round(regional_rsst_avg_2024, digits=2)) K")

#------------------------------------------------------
# functions to make plots: 
#
## plots only a region, defined by 'limits'
function fig_anom_reg_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:white,
        size=(600,300),
        )
    ax = Axis(f2[1,1]; #--> default plot is rectangular equidistant 
    #ax = GeoAxis(f2[1,1];
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
###
function fig_anom_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:white,
        size=(600,300),
        )
    ax = Axis(f2[1,1]; #--> default plot is rectangular equidistant 
    #ax = GeoAxis(f2[1,1];
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
function fig_3pan_plot(inp1,inp2,inp3,d1,d2,tit1,tit2,tit3, lon_min, lon_max, lat_min, lat_max)
    f2 = Figure(;
        #figure_padding=(5,5,10,10),
        backgroundcolor=:white,
        size=(300,1000),
        )
    ax1 = Axis(f2[1,1]; #--> default plot is rectangular equidistant 
        xticks = -180:20:180, 
        yticks = -90:10:90,
        ylabel="latitude",
        xticklabelsvisible = false,
        yticklabelsize = 20,
        ylabelsize = 20, #ylabelsize_bold = true,
        limits=(-120,0,0,40),
        title=tit1,
        titlesize=21.0,
        )
        bb = contourf!(ax1, d1, d2, inp1, 
             levels = range(-3, 3, length = 21), # rh
             colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        colsize!(f2.layout, 1, Aspect(1, 2.0))
        #ax1.xlabelfont = :bold
        #ax1.ylabelfont = :bold
        lines!(ax1, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        # Draw rectangle for the region of interest
        lines!(ax1, [lon_min, lon_max, lon_max, lon_min, lon_min], 
                   [lat_min, lat_min, lat_max, lat_max, lat_min],
               color = :black, linewidth = 3.0, label = "MDR")
        Colorbar(f2[1,2], bb)
    ax2 = Axis(f2[2,1]; #--> default plot is rectangular equidistant 
        xticks = -180:20:180, 
        yticks = -90:10:90,
        ylabel="latitude",
        yticklabelsize = 20,
        xticklabelsvisible = false,
        ylabelsize = 20, #ylabelsize_bold = true,
        limits=(-120,0,0,40),
        title=tit2,
        titlesize=21.0,
        )
        bb = contourf!(ax2, d1, d2, inp2, 
             levels = range(-3, 3, length = 21), # rh
             colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        #ax2.xlabelfont = :bold
        #ax2.ylabelfont = :bold
        lines!(ax2, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[2,2], bb)
    ax3 = Axis(f2[3,1]; #--> default plot is rectangular equidistant 
        xticks = -180:20:180, 
        yticks = -90:10:90,
        ylabel="latitude",
        xlabel="longitude",
        yticklabelsize = 20,
        xticklabelsize = 20,
        xlabelsize = 20,
        ylabelsize = 20, #ylabelsize_bold = true,
        limits=(-120,0,0,40),
        title=tit3,
        titlesize=21.0,
        )
        #ax3.xlabelfont = :bold
        #ax3.ylabelfont = :bold
        bb = contourf!(ax3, d1, d2, inp3, 
             levels = range(-3, 3, length = 21), # rh
             colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        lines!(ax3, GeoMakie.coastlines(), color = :black, linewidth=0.75)
    resize_to_layout!(f2)
    return f2
end
###

#tit="2005 ERA5 SST"
rsst1_m_rsst2 = rsst1_mn .- rsst2_mn
#fig = fig_anom_reg_plot(rsst1_m_rsst2[:,:],lon,lat,tit)
#fig1name=tag*"_rsst_reganom_2005C.png"
#save(fig1name, fig)
#
#tit="2013 ERA5 SST"
#fig = fig_anom_reg_plot(rsst2_mn[:,:],lon,lat,tit)
#fig1name=tag*"_rsst_reg_2013C.png"
#save(fig1name, fig)
#
rsst3_m_rsst2 = rsst3_mn .- rsst2_mn
#tit="2024 ERA5 SST"
#fig = fig_anom_reg_plot(rsst3_m_rsst2[:,:],lon,lat,tit)
##fig = fig_anom_reg_plot(rsst3_mn[:,:],lon,lat,tit)
#fig1name=tag*"_rsst_reganom_2024C.png"
#save(fig1name, fig)
#
data_anom  = data_2_plot .- data_2_plot2
data_anom2 = data_2_plot3 .- data_2_plot2
#tit="ERA5 SST anom 2005 - 2013"
#fig = fig_anom_plot(data_anom[:,:],lon,lat,tit) 
#fig1name=tag*"_sst_2005anomC.png"
#save(fig1name, fig)
#tit="ERA5 SST anom 2024 - 2013"
#fig = fig_anom_plot(data_anom2[:,:],lon,lat,tit) 
#fig1name=tag*"_sst_2024anomC.png"
#save(fig1name, fig)

tit1="RSST 2013"
tit2="Anomolous RSST 2005"
tit3="Anomolous RSST 2024"
fig = fig_3pan_plot(rsst2_mn[:,:],rsst1_m_rsst2[:,:],rsst3_m_rsst2[:,:], lon, lat, tit1, tit2, tit3, lon_min, lon_max, lat_min, lat_max)
#function fig_tot_plot_with_region(inpv, d1, d2, tit, lon_min, lon_max, lat_min, lat_max)
fig3panName=tag*"_3pan_blah"*res_tag*".png"
save(fig3panName, fig)
