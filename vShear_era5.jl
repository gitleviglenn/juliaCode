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
#using Plotly
using GeoMakie
using NCDatasets
using Statistics
using HypothesisTests
using TypedTables
using DataFrames

include("ensoFuncs.jl")

#
# create the indices that correspond to Nino and Nina years/months
# NH
#ninoyears = [18 54 90 150 174 234 306 402]
#ninayears = [102 114 210 246 318 366 378 390]
# SH
####ninoyears = [30 42 66 102 162 246 318 354]
ninoyears = [23 35 59 96 155 239 311 347]
ninayears = [107 119 215 251 263 335 371 383]
#####ninayears = [114 126 222 258 270 342 378 390]
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
timeENSO  = collect(1.:1:numfields);
era5_sst_1990th2023_360x180timeENSO  = collect(1.:1:48);
# span +/-40 degrees
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
VWS_high_full  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
VWS_high       = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])

endi = 48
# low should be an array that contains the timesteps representing the negative phase of ENSO
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
# when should i take the time average?    It doesn't seem to matter.

VWS_high_tmn = mean(VWS_high_full, dims=3)
VWS_low_tmn = mean(VWS_low_full, dims=3)

vws_ninoMnina = VWS_high_full .- VWS_low_full
vws_comp_mn   = mean(vws_ninoMnina, dims=3)
println("size of vsw_ninoMnina",size(vws_ninoMnina))

VWS_tot_tm = mean(VWS_tot, dims = 3)

data_2_plot_anom = VWS_high_tmn - VWS_low_tmn
data_2_plot_tot = VWS_tot_tm

#--------------------------------------------------------------------------------
# calculate significance:
X2 = Float64.(timeENSO);
d1 = 360
d2 = 81
pvalgrid  = Array{Union{Missing, Float64}, 2}(undef, d1, d2)
pvallow   = Array{Union{Missing, Float64}, 2}(undef, d1, d2)
pvalhigh  = Array{Union{Missing, Float64}, 2}(undef, d1, d2)
intVal    = Array{Union{Missing, Float64}, 2}(undef, d1, d2)

for i in 1:360 # dims[1] longitudes
    for j in 1:81 # dims[2] latitudes
        Yind = vws_ninoMnina[i,j,:];
        if any(ismissing, Yind)
            #println("missing value found ")
        else
        # need to somehow check if Yind contains missing data, and if so, then don't proceed to the lm step.  
            Ytemp = Float64.(Yind)
            tab   = Table(X = X2, Y = Ytemp);
            ols_temp = OneSampleTTest(Ytemp)
            pvalgrid[i,j] = pvalue(ols_temp)
            rang = confint(ols_temp, level = 0.95, tail = :both);
            pvallow[i,j]  = rang[1]
            pvalhigh[i,j] = rang[2]
            # test to make sure that 0 is not contained between the upper/lower 95% range
            intVal[i,j] = sign(pvalhigh[i,j]) + sign(pvallow[i,j])
        end
    end
end

arrN = reshape(pvalgrid, (d1*d2));
sortedP = sort(arrN)

lengthP = d1*d2
pAxis = collect(1:1:lengthP);
wilks0p05 = (0.05/lengthP) .* pAxis
wilks0p1 = (0.1/lengthP) .* pAxis
p0p05 = zeros(lengthP) .+ 0.05;
wilksArray = zeros(lengthP);
#---------------------------------
fig3 = Figure(;
    size = (400,400),
    )
ax = Axis(fig3[1,1];
    xlabel = "p-value rank i",
    ylabel = "p-value",
    title  = "Significance?",
    limits=(0,20000,0,0.1),
    )
#lines!(C,ts_roni_sm, linestyle = :solid)
lines!(ax,pAxis[:],sortedP[:], color = :black)
lines!(ax,pAxis[:],p0p05[:])
lines!(ax,pAxis[:],wilks0p05[:])
fig3
save("sig_vws_Wilks.png", fig3, px_per_unit=6.0)
#---------------------------------

for i in 1:20000
    wilksArray[i] = sortedP[i] - wilks0p05[i]
end

#psign = sign.(sortedP)
psign = sign.(wilksArray)
diffs = diff(psign)
psign_ind = findall(!iszero, diffs)
println("change of sign should be intersection point of sortedP and wilks0p05: ",psign_ind)
println("pvalue at intersection of sorted pvalues and Wilks FDR pvalue is: ",sortedP[psign_ind[1]])

# this is the p-value based on the FDR (false detection rate) described in Wilkes, 2016
pFDR = sortedP[psign_ind[1]]

# create a boolean array of points that indicate statistical significance
woman      = falses(d1,d2);
c1 = coalesce.(pvalgrid, false);
c2 = coalesce.(intVal, false);
for i in 1:d1
    for j in 1:d2
        #woman[i,j] = c1[i,j] < 0.05 && c2[i,j] != 0.0
        woman[i,j] = c1[i,j] < pFDR && c2[i,j] != 0.0
    end
end
# grab all the points that are labelled as 'true' or 1.  
points = findall(x -> x == 1, woman);
# points is now a CartesianIndex object, which cannot be directly used with scatter.  
# we have to create coordinate arrays from points: 
x_coords = [idx.I[1] for idx in points];
y_coords = [idx.I[2] for idx in points];

# shift index arrays to correspond to trad lat/lon defs
x_lats = x_coords .- 180;
#x_lats = x_coords
#y_lons = y_coords .- 90; # shifting the south point from 0 to -90
y_lons = y_coords .- 40; # shifting the south point from 0 to -40

println("=========================")
#println("x_lats are: ",x_lats[:])
println("=========================")
#println("y_lons are: ",y_lons[:])
println("=========================")
#println("lons are: ",lon[:])
println("=========================")
#println("lats are: ",lat[:])
println("=========================")

#--------------------------------------------------------------------------------

##

#--------------------------------------------
f4 = Figure(;
    figure_padding=(5,5,10,10),
    backgroundcolor=:white,
    size=(900,400),
    )
ax = GeoAxis(f4[1,1];
    xticks = -180:30:180, 
    #xticks = 0:30:360, 
    yticks = -90:30:90,
    ylabel="latitude",
    xlabel="longitude",
    limits=(-180,180,-40,40),
    title="VWS[m/s] composite El Nino - La Nina",
    xticklabelsize = 22, # 14,16 are pretty reasonable sizes
    yticklabelsize = 22, # 22 used for 8 panel figure that needs larger font
    )
    bb = contourf!(ax, lon, lat, vws_comp_mn[:,:,1], 
         #levels = range(0, 50, length = 25), # tos
         levels = range(-10, 10, length = 21), # rh
         #colormap = :Blues_8,
         #colormap = :broc,
         colormap = :bam,
         #colormap = :batlow,
         #colormap = :vik,
         extendlow = :auto, extendhigh = :auto
    )
    lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
    Colorbar(f4[1,2], bb)
scatter!(ax, x_lats, y_lons, marker = :circle, markersize=2, color = :black)
f4

#vws_comp_mn
#fig = fig_anom_plot(vws_comp_mn[:,:],lon,lat,tit,levs)
save("vws_fig_SH_Wilks.png", f4, px_per_unit=6.0)


