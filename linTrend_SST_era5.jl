#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# linTrend_SST_era5.jl
#    --> modified from linTrendSST_OctNov_era5.jl
# plot regional map of linear trend of SST
#
# a plot of the statistical significance according to the methodology of Wilks, 2016 is also created.
# the default behavior is to use a p value of 0.05.  0.1 can also be used if desired.   look for the 
# code labelled as: Wilks threshold
#
# label grid points that have trends with statistical significance at the 0.05 level following Wilks,
# 2016.
#
# when changing the data used with this script, think carefully about the units.  currently the 
# data is scaled by a factor of 10 to give C/decade.   this is because each time step is one year. 
# if the data is monthly then a scaling factor of 120 is necessary.   for october november data, 
# which is monthly, but for which each 'year' only contains 2 months, the scaling factor for units
# of C/decade is 20.  
# 
# currently it is assumed that each time step is 1 degree of freedom.   If this is not the case then 
# look at the code block labelled 'DOF' and adjust accordingly.
#
# levi silvers                                                                  dec 2025
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics
using GLM
using TypedTables
using Polynomials
using ShiftedArrays

include("ensoFuncs.jl")

#colormap = :batlow,
#colormap = :bam, # default for shear plot (greens and pinks)
#colormap = :seismic, # colors are a bit harsh
#colormap = :vik, # default for redish bluish for relative SST
#colormap = :BrBg, # better for RH  browns and greens
#colormap = :roma,

function fig_plot_stats(inpv,d1,d2,tit,levs,cbar,xPts,yPts)
    funfig = Figure(;
        fontsize = 20,
        figure_padding=(5,5,10,10),
        #figure_padding=(10,15,10,10),
        backgroundcolor=:white,
        size=(900,400),
        #size=(600,300), # this increases tickfont size, but doesn't print -30S!!#$%
        )
    ax = GeoAxis(funfig[1,1];
        dest = "+proj=cea",
        xticks = -180:20:180,
        #xticks = -160:20:180,
        #xticks = 0:30:360, 
        yticks = -90:20:90,
        #yticks = -60:20:60,
        #yticks = -60:20:60,
        xlabel="longitude",
        ylabel="latitude",
        #limits=(-160,20,-20,40),
        limits=(-180,180,-60,60),
        title=tit,
        xticklabelsize = 12, # 14,16 are pretty reasonable sizes
        yticklabelsize = 15, # 22 used for 8 panel figure that needs larger font
        )
        bb = contourf!(ax, d1, d2, inpv,
          levels = levs,
          colormap = cbar,
          extendlow = :auto, extendhigh = :auto
        )
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(funfig[1,2], bb)
        scatter!(ax, xPts, yPts, marker = :utriangle, markersize=2, color = :black)
    return funfig
end

path="/Users/C823281551/data/ERA5/"

# this file contains ERA5 data, two months per year.
#filein1  = path*"era5_sst_1979th2024_OctNov_360x180.nc"

# this file contains ERA5 data, the yearly mean has already been computed with NCO
# so each timestep corresponds to 1 year.   
filein1  = path*"era5_sst_annmn_1990th2024_360x180.nc"

tag = "ERA5"
data  = NCDataset(filein1)

lat = data["lat"]
lon = data["lon"]
tme = data["valid_time"]
 
dim_var  = data["sst"]

timelen = 35 
#timelen = 46 
##timeAxis1 = collect(1:1:92);
#timeAxis  = collect(1:1:46);
timeAxis  = collect(1:1:timelen);

sst_var = data["sst"];
dims    = size(dim_var)

print("size of sst Array is: ",size(dim_var))

# -----------------------------------------
# degrees of freedom 

# average the months for each year to result in arrays with 8 time stamps for 8 years.
# this allows for statistically testing with 7 degrees of freedome instead of 47.

sst_yrs    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], timelen)

# the code below is very important for specifying the degrees of freedom.  for example
# with the case of looking at variables during October and November, we assumed that 
# each year was independent, but not each month.   So we averaged over the two months
# of each year
#
# DOF : degrees of freedom
#jj = 1
#for i in 1:2:dims[3]-1 # assume each year is an independant degree of freedom
##for i in 1:2:timelen-1 # assume each year is an independant degree of freedom
#    test1 = sst_var[:,:,i:i+1];
#    sst_yrs[:,:,jj]  = mean(test1, dims=3)
#    global jj += 1
#end

sst_yrs = sst_var

println("dims[1] and dims[2] are: ",dims[1]," and ",dims[2])
println("the size of sst_yrs is: ",size(sst_yrs))
println("size of timeAxis is: ",size(timeAxis))

# -----------------------------------------

# these arrays store the linear regression values at each grid point
bgrid   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])

#sst_trop = Array{Union{Missing, Float64}, 2}(undef, dims[1], 41)
sst_trop_mn = Array{Union{Missing, Float64}, 1}(undef, 35)

#    sst_yrs[:,:,jj]  = mean(test1, dims=3)
for i in 1:dims[1]
    for j in 1:dims[2]
        agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis[:],sst_yrs[i,j,:])
    end
end

#for t in 1:35
# sst_trop_mn[t] = mean(skipmissing(agrid[:,:,t]))
#end 
sst_trend_mn = mean(skipmissing(agrid))

println("tropical mean SST trends are: ",sst_trend_mn*10)

#sst_mn_v1 = mean(sst_yrs, dims = 1)
#println("size of sst_mn_v1 is: ",size(sst_mn_v1))
#sst_mn_v2 = mean(sst_mn_v1, dims = 2)
#println("size of sst_mn_v2 is: ",size(sst_mn_v2))
#sst_mn_v3 = mean(sst_mn_v2, dims = 3)
#println("size of sst_mn_v3 is: ",size(sst_mn_v3))

#----------------------------------------------------------------------------
# work on significance
#----------------------------------------------------------------------------
#d1 = 360
d1 = dims[1]
d2 = 81
pvalgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pvallow   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pvalhigh  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
intVal    = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])

println("sst_var values at 1 grid point are: ",sst_var[10,90,:])

# why is the lm model needed here, but the OneSampleTTest needed for the 
# composite figures?  

for i in 1:360 # dims[1] longitudes
  for j in 1:180 # dims[2] latitudes
  #for j in 50:130 # dims[2] latitudes
    X2 = Float64.(timeAxis);
    #Yind = vws_yrs[i,j,:];
    Yind = sst_yrs[i,j,:];
    #Yind = VWS_tot[i,j,:];
    #Yind = rh_lev[i,j,:];
    #Yind = sst_var[i,j,:];
    if any(ismissing, Yind)
      #println("missing value found ")
    else
      #println("missing value not found ")
      # need to somehow check if Yind contains missing data, and if so, then don't proceed to the lm step.  
      Ytemp = Float64.(Yind)
      #ols_temp = OneSampleTTest(Ytemp)
      #pvalgrid[i,j-49] = pvalue(ols_temp)
      #rang = confint(ols_temp, level = 0.95, tail = :both);
      #pvallow[i,j-49]  = rang[1]
      #pvalhigh[i,j-49] = rang[2]
      #intVal[i,j-49] = sign(pvalhigh[i,j-49]) + sign(pvallow[i,j-49])
      ##
      tab   = Table(X = X2, Y = Ytemp);
      ols_temp = lm(@formula(Y ~ X), tab)
      testVals = coeftable(ols_temp)
      #pvaltemp = testVals.cols[4][2]
      pvalgrid[i,j] = testVals.cols[4][2]
      pvallow[i,j]  = testVals.cols[5][2]
      pvalhigh[i,j] = testVals.cols[6][2]
      intVal[i,j] = sign(pvalhigh[i,j]) + sign(pvallow[i,j])
    end
  end
end

arrN = reshape(pvalgrid, (dims[1]*dims[2]));
sortedP = sort(arrN)

#lengthP = 81*360;
lengthP = 121*360;
fullength = 360*180

pAxis = collect(1:1:fullength);
wilks0p05 = (0.05/(lengthP)) .* pAxis
wilks0p1 = (0.1/(lengthP)) .* pAxis
p0p05 = zeros(fullength) .+ 0.05;
wilksArray = zeros(fullength);

#println("Wilks value for 0.05 is: ",wilks0p05)
#println("Wilks value for 0.1 is: ",wilks0p1)

#----------------------
xmax = 40000

fig3 = Figure(;
    size = (600,400),
    )
ax = Axis(fig3[1,1];
    xlabel = "p-value rank i",
    ylabel = "p-value",
    title  = "Significance?",
    limits=(0,xmax,0,0.1),
    )
lines!(ax,pAxis[:],sortedP[:], label = "sorted p-vales", color = :black)
lines!(ax,pAxis[:],p0p05[:], label = "p = 0.05")
lines!(ax,pAxis[:],wilks0p05[:], label = "Wilks 0.05")
lines!(ax,pAxis[:],wilks0p1[:], label = "Wilks 0.1")
fig3[1,2] = Legend(fig3, ax, " ", framevisible = false)
fig3
save("sig_Wilks_sst.png", fig3, px_per_unit=6.0)
#----------------------

println("size of sortedP is: ",size(sortedP))
println("size of wilks0p05is: ",size(wilks0p05))

# Wilks Threshold:
for i in 1:xmax
    wilksArray[i] = sortedP[i] - wilks0p05[i] # using a p value of 0.05
    #wilksArray[i] = sortedP[i] - wilks0p1[i]  # using a p value of 0.1
end

#psign = sign.(sortedP)
psign = sign.(wilksArray)
diffs = diff(psign)
psign_ind = findall(!iszero, diffs)
println("change of sign should be intersection point of sortedP and wilks0p05: ",psign_ind)
println("pvalue at intersection of sorted pvalues and Wilks FDR pvalue is: ",sortedP[psign_ind[1]])

# this is the p-value based on the FDR (false detection rate) described in Wilkes, 2016
pFDR = sortedP[psign_ind[1]]

# create a bitmatrix to identify statistically significane grid points
woman      = falses(dims[1],dims[2]);
c1 = coalesce.(pvalgrid, false);
c2 = coalesce.(intVal, false);
for i in 1:360
    for j in 1:180
        #woman[i,j] = pvalgrid[i,j] < 0.05 #&& intVal[i,j] != 0.0
        #woman[i,j] = c1[i,j] < 0.05 && c2[i,j] != 0.0
        woman[i,j] = c1[i,j] < pFDR && c2[i,j] != 0.0
    end
end

# grab all the points that are labelled as 'true' or 1.  
points = findall(x -> x == 1, woman); # points has dims of (24065,)
# points is now a CartesianIndex object, which cannot be directly used with scatter.  
# we have to create coordinate arrays from points: 
x_coords = [idx.I[1] for idx in points]; # x_coords has dims of (24065,)
y_coords = [idx.I[2] for idx in points];

# shift index arrays to correspond to trad lat/lon defs
xshift = 181
#yshift = 92
yshift = 91

x_lats_sst = x_coords .- xshift;
y_lons_sst = y_coords .- yshift;

#println("***********************")

# print figure panels

## sst levels
#levs = range(-0.5, 0.5, length = 21)
levs = range(-0.3, 0.3, length = 19)
#blah = fig_plot_stats(agrid.*20,lon,lat,"Annual SST linear trends (C)/decade",levs,:vik, x_lats_sst, y_lons_sst)
blah = fig_plot_stats(agrid.*10,lon,lat,"Annual SST linear trends (C)/decade",levs,:vik, x_lats_sst, y_lons_sst)
save("era5_SST_LinTrend_Test_AnnMn_1979th2024_35dof.png", blah, px_per_unit=6.0)
#
## VWS levels

## PI levels
















