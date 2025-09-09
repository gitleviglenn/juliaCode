#-----------------------------------------------------------------------------------------------
# rh_era5.jl
#
#-----------------------------------------------------------------------------------------------

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics
using GLM
using TypedTables

include("ensoFuncs.jl")
#
tag = "Dude"
# create the indices that correspond to Nino and Nina years/months
#fig2name = tag*"_rh_SH_test.png"
# NH
fig2name = tag*"_rh_nino_comp_NH_test.png"
ninoyears = [18 54 90 150 174 234 306 402]
ninayears = [102 114 210 246 318 366 378 390]
# SH
#fig2name = tag*"_rh_nino_comp_SH_test.png"
#ninoyears = [23 35 59 96 155 239 311 347]
#ninayears = [107 119 215 251 263 335 371 383]
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

filein  = path*"era5_rh_1990th2023_360x80.nc"
fileclimo = path*"era5_rh_1990th2023_360x80_climo_climatology.nc"

tag = "ERA5"
data   = NCDataset(filein)

numfields = 48 
numall    = 408
timeENSO  = collect(1.:1:numfields);

lat = data["lat"]
lon = data["lon"]
lev = data["pressure_level"]

rh_var_full  = data["r"]
rh_var_climo = data["r"]
dims = size(rh_var_full)

# remove climatology from rh_var?   Use something like this: 
#  # computed ts for nino 3.4 minus seasonal cycle
#  global ts_rmn_nsc = zeros(tlength) # time series of nino3.4 without seasonal cycle
#  for i in 1:12:tlength
#    ts_rmn_nsc[i]    = nino34_ts[i]    - ss[1]
#    ts_rmn_nsc[i+1]  = nino34_ts[i+1]  - ss[2]
#    ts_rmn_nsc[i+2]  = nino34_ts[i+2]  - ss[3]
#    ts_rmn_nsc[i+3]  = nino34_ts[i+3]  - ss[4]
#    ts_rmn_nsc[i+4]  = nino34_ts[i+4]  - ss[5]
#    ts_rmn_nsc[i+5]  = nino34_ts[i+5]  - ss[6]
#    ts_rmn_nsc[i+6]  = nino34_ts[i+6]  - ss[7]
#    ts_rmn_nsc[i+7]  = nino34_ts[i+7]  - ss[8]
#    ts_rmn_nsc[i+8]  = nino34_ts[i+8]  - ss[9]
#    ts_rmn_nsc[i+9]  = nino34_ts[i+9]  - ss[10]
#    ts_rmn_nsc[i+10] = nino34_ts[i+10] - ss[11]
#    ts_rmn_nsc[i+11] = nino34_ts[i+11] - ss[12]
#  end

rh_var  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numall)

rh_var = rh_var_full;

for i in 1:360 # dims[1] longitudes
  for j in 1:81 # dims[2] latitudes
    for k in 1:12:numall
      #rh_var[i,j,k] = rh_var_full[i,j,k] .- rh_var_climo[i,j,1]
      #rh_var[i,j,k] = rh_var_full[i,j,k+1] - rh_var_climo[i,j,2]
      #rh_var[i,j,k] = rh_var_full[i,j,k+2] - rh_var_climo[i,j,3]
      #rh_var[i,j,k] = rh_var_full[i,j,k+3] - rh_var_climo[i,j,4]
      #rh_var[i,j,k] = rh_var_full[i,j,k+4] - rh_var_climo[i,j,5]
      #rh_var[i,j,k] = rh_var_full[i,j,k+5] - rh_var_climo[i,j,6]
      #rh_var[i,j,k] = rh_var_full[i,j,k+6] - rh_var_climo[i,j,7]
      #rh_var[i,j,k] = rh_var_full[i,j,k+7] - rh_var_climo[i,j,8]
      #rh_var[i,j,k] = rh_var_full[i,j,k+8] - rh_var_climo[i,j,9]
      #rh_var[i,j,k] = rh_var_full[i,j,k+9] - rh_var_climo[i,j,10]
      #rh_var[i,j,k] = rh_var_full[i,j,k+10] - rh_var_climo[i,j,11]
      #rh_var[i,j,k] = rh_var_full[i,j,k+11] - rh_var_climo[i,j,12]
    end
  end
end

lat1      = 1
lat2      = 81

rh_1            = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
rh_2            = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
rh_low          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
rh_high         = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
rh_tot          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numall)

endi = 48
# low is an array that contains the timesteps representing the negative phase of ENSO
# high is an array that contains the timesteps representing the positive phase of ENSO
low = ninaInd
high = ninoInd

# select the case determining if rh will be plotted at one level or at the average of two
pcase = 0# pcase = 0 corresponds to plotting rh on 1 level
# calculate total rh field, for all times
level=2 # level 2 should correspond to the 700 hPa pressure level. 
if pcase == 0
  for i in 1:408
    rh_tot[:,lat1:lat2,i]   = rh_var[:,lat1:lat2,level,i] 
  end
  for i in 1:endi
    # low values, la nina
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,level,low[i]]
    rh_low[:,lat1:lat2,i]   = rh_2[:,lat1:lat2,1,i] 
    # high values, el nino
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,level,high[i]]
    rh_high[:,lat1:lat2,i]  = rh_2[:,lat1:lat2,1,i] 
  end
else
  for i in 1:408
    rh_tot[:,lat1:lat2,i]   = (rh_var[:,lat1:lat2,2,i] + rh_var[:,lat1:lat2,1,i])./2
  end
  for i in 1:endi
    # low values, la nina
    rh_1[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,1,low[i]]
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,2,low[i]]
    rh_low[:,lat1:lat2,i]   = (rh_2[:,lat1:lat2,1,i] + rh_1[:,lat1:lat2,1,i])./2
    # high values, el nino
    rh_1[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,1,high[i]]
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,2,high[i]]
    rh_high[:,lat1:lat2,i]  = (rh_2[:,lat1:lat2,1,i] + rh_1[:,lat1:lat2,1,i])./2
  end
end 


rh_tot_tmn  = mean(rh_tot, dims=3)
rh_high_tmn = mean(rh_high, dims=3)
rh_low_tmn  = mean(rh_low, dims=3)

rh_ninoMnina = rh_high .- rh_low
rh_comp_mn   = mean(rh_ninoMnina, dims=3)

data_2_plot_tot  = rh_tot_tmn
data_2_plot_anom = rh_high_tmn - rh_low_tmn

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
        Yind = rh_ninoMnina[i,j,:];
        if any(ismissing, Yind)
            #println("missing value found ")
        else
        # need to somehow check if Yind contains missing data, and if so, then don't proceed to the lm step.  
            Ytemp = Float64.(Yind)
            tab   = Table(X = X2, Y = Ytemp);
            ols_temp = lm(@formula(Y ~ X), tab)
            testVals = coeftable(ols_temp)
            pvalgrid[i,j] = testVals.cols[4][2]
            pvallow[i,j]  = testVals.cols[5][2]
            pvalhigh[i,j] = testVals.cols[6][2]
            # test to make sure that 0 is not contained between the upper/lower 95% range
            intVal[i,j] = sign(pvalhigh[i,j]) + sign(pvallow[i,j])
        end
    end
end

arrN = reshape(pvalgrid, (360*81));
sortedP = sort(arrN)

lengthP = 360*81
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
    limits=(0,10000,0,0.1),
    )
#lines!(C,ts_roni_sm, linestyle = :solid)
lines!(ax,pAxis[:],sortedP[:], color = :black)
lines!(ax,pAxis[:],p0p05[:])
lines!(ax,pAxis[:],wilks0p05[:])
fig3
save("sig_vws_Wilks.png", fig3, px_per_unit=6.0)
#---------------------------------

# create boolean array to plot markers where grid points are significant: 
for i in 1:10000
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

#---------------------------------
# below is code for original figure: 
tit="ERA5 RH Composite (%)"
levs = range(-15., 15., length = 21)
#fig2name = tag*"_rh_nino_comp_SHmn.png"
#fig = fig_tot_plot(data_2_plot_tot[:,:,1],lon,lat,tit)
fig = fig_anom_plot(data_2_plot_anom[:,:,1],lon,lat,tit,levs)
#lonlon=[1:360;]

#fig = fig_anom_plot(data_2_plot_anom[:,:],lon,lat,tit)
save(fig2name, fig)

#-----------------------------------------------------------------------------------------------------------
f4 = Figure(;
    figure_padding=(5,5,10,10),
    backgroundcolor=:white,
    size=(900,400),
    )
ax = GeoAxis(f4[1,1];
    xticks = -180:30:180, 
    #xticks = 0:30:360, 
    yticks = -40:20:40,
    ylabel="latitude",
    xlabel="longitude",
    limits=(-180,180,-40,40),
    title=tit,
    xticklabelsize = 16, # 14,16 are pretty reasonable sizes
    yticklabelsize = 16, # 22 used for 8 panel figure that needs larger font
    )
    bb = contourf!(ax, lon, lat, rh_comp_mn[:,:,1], 
         #levels = range(0, 50, length = 25), # tos
         levels = range(-15, 15, length = 21), # rh
         #colormap = :Blues_8,
         #colormap = :broc,
         #colormap = :bam,
         colormap = :BrBg,
         #colormap = :batlow,
         #colormap = :vik,
         extendlow = :auto, extendhigh = :auto
    )
    lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
    Colorbar(f4[1,2], bb)
scatter!(ax, x_lats, y_lons, marker = :utriangle, markersize=4, color = :black)
f4

#vws_comp_mn
#fig = fig_anom_plot(vws_comp_mn[:,:],lon,lat,tit,levs)
save("rh_fig.png", f4, px_per_unit=6.0)

