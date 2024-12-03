#-----------------------------------------------------------------------------------------------
# vShearMap.jl
#-----------------------------------------------------------------------------------------------
using CairoMakie
using GeoMakie
using NCDatasets

path="/Users/C823281551/data/"
modelp="MPI-ESM1-2-LR" 
#filehur  = path*"cmip6/MPIESM/hur_Amon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216.nc" 

#filein  = path*"cmip6/MPIESM/hur_Amon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_regridded3.nc" 
#filein  = path*"cmip6/MPIESM/tos_Omon_"*modelp*"_ssp585_r1i1p1f1_gn_20150116-21001216_pm40b.nc" 

filein   = path*"cmip6/CESM2/ua_Amon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215_360x180.nc"
filein2  = path*"cmip6/CESM2/va_Amon_CESM2_ssp585_r4i1p1f1_gn_20150115-21001215_360x180.nc"

data   = NCDataset(filein)
data2  = NCDataset(filein2)

data.attrib

lat  = data["lat"]
lon  = data["lon"]
lev  = data["plev"]
tme = data["time"]

u_var  = data["ua"] # hur(time, plev, lat, lon)
u_var  = data["ua"] # hur(time, plev, lat, lon)
v_var  = data2["va"] # hur(time, plev, lat, lon)
v_var  = data2["va"] # hur(time, plev, lat, lon)
tit = "wind: CESM2"

# tropics in CESM2
lat1=50
lat2=130
plev=2

tim = 506 # time step at which data will be plotted.

data_hi  = u_var[:,lat1:lat2,2,tim]
data_low = u_var[:,lat1:lat2,1,tim]
uSh      = data_hi - data_low
data_hi  = v_var[:,lat1:lat2,2,tim]
data_low = v_var[:,lat1:lat2,1,tim]
vSh      = data_hi - data_low

#data_2_plot = data_hi - data_low
data_2_plot = sqrt.(uSh.^2 .+ vSh.^2)
#data_2_plot = sqrt.(uSh.^2)

function fig_1_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(600,300),
        )
    ax = Axis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        limits=(-180,180,-40,40),
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             levels = range(0, 50, length = 25), # tos
             #levels = range(0, 100, length = 10), # rh
             #colormap = :Blues_8,
             #colormap = :navia,
             colormap = :batlow,
             #colormap = :vik,
        )
        lines!(ax, GeoMakie.coastlines())
        Colorbar(f2[1,2], bb)
    return f2
end

fig = fig_1_plot(data_2_plot,lon,lat[lat1:lat2],tit)


