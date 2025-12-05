#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# linTrend_SST_era5.jl
#    --> modified from linTrendSST_OctNov_era5.jl
# plot regional map of linear trend of SST
#
# label grid points that have trends with statistical significance at the 0.05 level following Wilks,
# 2016.
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
        #xticks = -180:20:180,
        xticks = -170:20:170,
        #xticks = 0:30:360, 
        #yticks = -90:30:90,
        yticks = -60:20:60,
        xlabel="longitude",
        ylabel="latitude",
        limits=(-160,20,-20,40),
        title=tit,
        xticklabelsize = 16, # 14,16 are pretty reasonable sizes
        yticklabelsize = 16, # 22 used for 8 panel figure that needs larger font
        )
        bb = contourf!(ax, d1, d2, inpv,
          levels = levs,
          colormap = cbar,
          extendlow = :auto, extendhigh = :auto
        )
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(funfig[1,2], bb)
        scatter!(ax, xPts, yPts, marker = :utriangle, markersize=4, color = :black)
    return funfig
end

path="/Users/C823281551/data/ERA5/"

filein1  = path*"era5_sst_1979th2024_OctNov_360x180.nc"

tag = "ERA5"
data  = NCDataset(filein1)

lat = data["lat"]
lon = data["lon"]
tme = data["valid_time"]
 
dim_var  = data["sst"]

timeAxis1 = collect(1:1:92);
timeAxis  = collect(1:1:46);

sst_var = data["sst"];
dims    = size(dim_var)

print("size of sst Array is: ",size(dim_var))

# -----------------------------------------
# degrees of freedom 

# average the months for each year to result in arrays with 8 time stamps for 8 years.
# this allows for statistically testing with 7 degrees of freedome instead of 47.

sst_yrs    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 46)

jj = 1
#for i in 1:6:48-5
for i in 1:2:dims[3]-1 # assume each year is an independant degree of freedom
    test1 = sst_var[:,:,i:i+1];
    sst_yrs[:,:,jj]  = mean(test1, dims=3)
    global jj += 1
end
# -----------------------------------------

# these arrays store the linear regression values at each grid point
bgrid   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])

for i in 1:dims[1]
    for j in 1:dims[2]
        agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis[:],sst_yrs[i,j,:])
    end
end

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
      println("missing value found ")
    else
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

lengthP = 81*360;
fullength = 360*180

pAxis = collect(1:1:fullength);
wilks0p05 = (0.05/(lengthP)) .* pAxis
wilks0p1 = (0.1/(lengthP)) .* pAxis
p0p05 = zeros(fullength) .+ 0.05;
wilksArray = zeros(fullength);

println("Wilks value for 0.05 is: ",wilks0p05)
println("Wilks value for 0.1 is: ",wilks0p1)

#----------------------
fig3 = Figure(;
    size = (400,400),
    )
ax = Axis(fig3[1,1];
    xlabel = "p-value rank i",
    ylabel = "p-value",
    title  = "Significance?",
    limits=(0,20000,0,0.1),
    )
lines!(ax,pAxis[:],sortedP[:], color = :black)
lines!(ax,pAxis[:],p0p05[:])
lines!(ax,pAxis[:],wilks0p05[:])
fig3
save("sig_Wilks_sst.png", fig3, px_per_unit=6.0)
#----------------------

for i in 1:20000
    wilksArray[i] = sortedP[i] - wilks0p05[i]
    #wilksArray[i] = sortedP[i] - wilks0p1[i]
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
levs = range(-0.5, 0.5, length = 21)
blah = fig_plot_stats(agrid.*20,lon,lat,"SST linear trends October-November ( C)/decade",levs,:vik, x_lats_sst, y_lons_sst)
save("era5_SST_LinTrend_Test_p0p05_1979th2024_region_46dof.png", blah, px_per_unit=6.0)
#
## VWS levels

## PI levels
















