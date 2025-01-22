#-----------------------------------------------------------------------------------------------
# vShear_era5.jl
#
# open era5 data for the u and v components of the wind.  compute the vertical wind shear (VWS)
#
# plot: mean VWS over entire period, mean VWS for El Nino, and mean VWS for La Nina 
#
# ERA5 data was downloaded using this script: 
# getERA5_uv_1990th2023.py
#
# ERA5 data was remapped like this: 
# cdo -L remapbil,mygridTropics era5_uv_1990th2023.nc era5_uv_1990th2023_360x80.nc
#
# NH:       Calculate for June through November
# Nino:     1991, 1994, 1997, 2002, 2004, 2009, 2015, 2023
# Neutral:  1990, 1992, 1993, 1995, 1996, 2000, 2001, 2003, 2005, 2006, 2008,
#           2011, 2012, 2013, 2014, 2017, 2018, 2019
# Nina:     1998, 1999, 2007, 2010, 2016, 2020, 2021, 2022
#
# SH:       Calculate for November through April
# Nino:     1992, 1993, 1995, 1998, 2003, 2010, 2016, 2019
# Neutral:  1990, 1991, 1994, 1996, 1997, 2001, 2002, 2004, 2005, 2007, 2009,
#           2012, 2013, 2014, 2015, 2017, 2020, 2023
# Nina:     1999, 2000, 2008, 2011, 2012, 2018, 2021, 2022
#
# 1990                1995                 2000                     2005                    2010
#   1  13  25  37  49  61  73  85  97  109  121  133  145 157  169  181  193  205  217  229  241 
#                      2015                    2020          
#  253  265  277  289  301  313  325  337  349  361  373  385  397
#
# levi silvers                                                                 jan 2025
#-----------------------------------------------------------------------------------------------
#
# create Nino indices for NH
# june through november of years 1991, 1994, 1997, 2002, 2004, 2009, 2015, and 2023
#
# create Nina indices for NH
# june through november of years 1998, 1999, 2007, 2010, 2016, 2020, 2021, and 2022
#
# assuming that i = 1 is january 1990, i = 13, janary 1991, i = 6, june 1990, i = 18, june 1991
# so i=18:23 would correspond to the TC season of 1991 during an el Nino
#
# si = 18
# i = si:1:si+5
# NH: enso +
# si = 18, 49+5, 85+5, 145+5, 169+5, 229+5, 301+5, 397+5 # nino starting indices
# NH: enso -
# si = 102 114 210 246 318 366 378 390
#
# --> these are wrong, because they need to be shifted than more than
# 5.   I think it should actually be 11, or maybe minus 1 or 2? 
# SH: enso +
# si = 30 42 66 102 162 246 318 354
# SH: enso -
# si = 114 126 222 258 270 342 378 390
# si = 109+5, 121+5, 217+5, 253+5, 265+5, 337+5, 373+5, 385+5
#
# o1 = 18:23
# o2 = 54:59
# o3 = 90:95
# o4 = 150:155
# o5 = 174:179
# o6 = 234:239
# o7 = 306:311
# o8 = 402:407
#
## create indices: input 8 years, output array/matric with 8*6 elements
#function create_indices(years)
#  for i in 1:8
#    si = years[i]
#    ensoInd=si:si+5
#  end
#  #o1=years[1]:years[1]+5
#  #ensoInd = [o1 o2 o3 o4 o5 o6 o7 o8]
#  return ensoInd
#end 
##
# si = 97+5, 109+5, 205+5, 241+5, 313+5, 361+5, 373+5, 385+5 # nina starting indices
#
# blue = 1:12:408
# julia> [i for i in blue] --> gives the list of january indices
# 1990              1995                 2000                    2005                     2010
# 1  13  25  37  49  61  73  85  97  109  121  133  145 157  169  181  193  205  217  229  241 
#                     2015                     2020          
#  253  265  277  289  301  313  325  337  349  361  373  385  397
#
#-----------------------------------------------------------------------------------------------

using CairoMakie
using GeoMakie
using NCDatasets

#
# create the indices that correspond to Nino and Nina years/months
# NH
#ninoyears = [18 54 90 150 174 234 306 402]
#ninayears = [102 114 210 246 318 366 378 390]
# SH
#ninoyears = [30 42 66 102 162 246 318 354]
ninoyears = [23 35 59 96 155 239 311 347]
#ninayears = [114 126 222 258 270 342 378 390]
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
#

path="/Users/C823281551/data/ERA5/"

filein  = path*"era5_uv_1990th2023_360x80.nc"
tag = "ERA5"
data   = NCDataset(filein)

lat = data["lat"]
lon = data["lon"]
lev = data["pressure_level"]
tme = data["valid_time"]
 
u_var = data["u"]
v_var = data["v"]

numfields = 48
numall    = 408
# probably will be spanning +/-40 degrees
lat1=1
lat2=81

dims = size(u_var)
u_1            = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
u_2            = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
uSh_tmp        = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
u_tot          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numall)
v_tot          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numall)
VWS_tot        = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numall)
v_1            = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
v_2            = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
vSh_tmp        = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
VWS_low_full   = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
#VWS_low_a      = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
#VWS_high_a     = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
VWS_high_full  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
VWS_high       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])

endi = 48
# low should be an array that contains the timesteps representing the negative phase of ENSO
#low = 1:1:400
#high = 1:1:400
low = ninaInd
high = ninoInd
# for low values of ENSO timeseries.   Also compute for high values.  
for i in 1:408
  # calculate the shear for all times
  u_tot[:,lat1:lat2,i]   = u_var[:,lat1:lat2,2,i] - u_var[:,lat1:lat2,1,i]
  v_tot[:,lat1:lat2,i]   = v_var[:,lat1:lat2,2,i] - v_var[:,lat1:lat2,1,i]
  VWS_tot[:,lat1:lat2,i] = sqrt.(u_tot[:,lat1:lat2,i].^2 .+ v_tot[:,lat1:lat2,i].^2)
end
for i in 1:endi
  # get low values
  u_1[:,lat1:lat2,1,i]   = u_var[:,lat1:lat2,1,low[i]]
  u_2[:,lat1:lat2,1,i]   = u_var[:,lat1:lat2,2,low[i]]
  uSh_tmp[:,lat1:lat2,i] = u_2[:,lat1:lat2,1,i] - u_1[:,lat1:lat2,1,i]
  v_1[:,lat1:lat2,1,i]   = v_var[:,lat1:lat2,1,low[i]]
  v_2[:,lat1:lat2,1,i]   = v_var[:,lat1:lat2,2,low[i]]
  vSh_tmp[:,lat1:lat2,i] = v_2[:,lat1:lat2,1,i] - v_1[:,lat1:lat2,1,i]
  VWS_low_full[:,lat1:lat2,i] = sqrt.(uSh_tmp[:,lat1:lat2,i].^2 .+ vSh_tmp[:,lat1:lat2,i].^2)
  # get high values
  u_1[:,lat1:lat2,1,i]   = u_var[:,lat1:lat2,1,high[i]]
  u_2[:,lat1:lat2,1,i]   = u_var[:,lat1:lat2,2,high[i]]
  uSh_tmp[:,lat1:lat2,i] = u_2[:,lat1:lat2,1,i] - u_1[:,lat1:lat2,1,i]
  v_1[:,lat1:lat2,1,i]   = v_var[:,lat1:lat2,1,high[i]]
  v_2[:,lat1:lat2,1,i]   = v_var[:,lat1:lat2,2,high[i]]
  vSh_tmp[:,lat1:lat2,i] = v_2[:,lat1:lat2,1,i] - v_1[:,lat1:lat2,1,i]
  VWS_high_full[:,lat1:lat2,i] = sqrt.(uSh_tmp[:,lat1:lat2,i].^2 .+ vSh_tmp[:,lat1:lat2,i].^2)
end
# when should i take the time average? 
#VWS_tot[:,lat1:lat2]  = sqrt.(u_tot[:,lat1:lat2,i].^2 .+ v_tot[:,lat1:lat2,i].^2)
#
#VWS_low_a[:,lat1:lat2]  = sqrt.(uSh_tmp[:,lat1:lat2,2].^2 .+ vSh_tmp[:,lat1:lat2,2].^2)
#VWS_high_a[:,lat1:lat2] = sqrt.(uSh_tmp[:,lat1:lat2,5].^2 .+ vSh_tmp[:,lat1:lat2,5].^2)

VWS_high_tmn = mean(VWS_high_full, dims=3)
VWS_low_tmn = mean(VWS_low_full, dims=3)

#VWS_total[:,lat1:lat2]  = sqrt.(u_2[:,lat1:lat2,1,i].^2 .+ vSh_tmp[:,lat1:lat2,2].^2)

#VWS_high_a[:,lat1:lat2] = sqrt.(uSh_tmp[:,lat1:lat2,5].^2 .+ vSh_tmp[:,lat1:lat2,5].^2)
#nino34_ts_mn = mean(filter(!isnan, skipmissing(nino34_full)))

#VWS_low  = mean(VWS_low_full, dims =3)
#VWS_high = mean(VWS_high_full, dims =3)

VWS_tot_tm = mean(VWS_tot, dims = 3)

#VWS_high = mean(filter(!isnan, skipmissing(VWS_high_full)), dims =3)
#VWS_high = mean(skipmissing(VWS_high_full), dims =3)
#tropmn_ts[i]=mean(skipmissing(trop_full[:,:,i]))-tropmn_ts_mn
#data_2_plot = VWS_high - VWS_low
#
#data_2_plot_anom = VWS_high_a - VWS_low_a

data_2_plot_anom = VWS_high_tmn - VWS_low_tmn
data_2_plot_tot = VWS_tot_tm

#
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
             levels = range(-20, 20, length = 21), # rh
             #colormap = :Blues_8,
             #colormap = :broc,
             colormap = :bam,
             #colormap = :batlow,
             #colormap = :vik,
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
             levels = range(0, 50, length = 20), # rh
             #colormap = :Blues_8,
             #colormap = :broc,
             colormap = :bam,
             #colormap = :batlow,
             #colormap = :vik,
             extendlow = :auto, extendhigh = :auto
        )
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[1,2], bb)
    return f2
end
tit="dummy"
fig2name = tag*"_vShear_nino_comp_SH.png"
#fig = fig_tot_plot(data_2_plot_tot[:,:,1],lon,lat,tit)
#fig = fig_anom_plot(data_2_plot_anom[:,:,1],lon,lat,tit)
fig = fig_anom_plot(data_2_plot_anom[:,:],lon,lat,tit)
save(fig2name, fig)
