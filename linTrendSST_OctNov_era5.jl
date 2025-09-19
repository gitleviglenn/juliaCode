# plot linear sst trend from October and November over 1979-2024
#
# levi silvers

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics
#using HypothesisTests
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
#filein1b = path*"era5_sst_1979th2024_OctNov.nc"
filein   = path*"MPI_ERA5_OctNov_full_output.nc"
filein2  = path*"era5_hur_OctNov_1979th2024_360x180.nc"
filein3  = path*"era5_uWind_OctNov_1979th2024_360x180.nc"
filein4  = path*"era5_vWind_OctNov_1979th2024_360x180.nc"

tag = "ERA5"
data1  = NCDataset(filein1)
#data1b = NCDataset(filein1b)
data2  = NCDataset(filein2)
data3  = NCDataset(filein3)
data4  = NCDataset(filein4)
data   = NCDataset(filein)

#lat = data["latitude"]
#lon = data["longitude"]
lat = data["lat"]
lon = data["lon"]
plev = data2["pressure_level"]
tme = data["valid_time"]
#latb = data1b["latitude"]
#lonb = data1b["longitude"]
 
dim_var  = data1["sst"]
#dim_varb = data1b["sst"]

timeAxis1 = collect(1:1:92);
timeAxis  = collect(1:1:46);

sst_var = data["sst"];
#sst_var_b = data1b["sst"];
pi_var  = data["vmax"];
dims    = size(dim_var)
#dim_high= size(dim_varb)
rh_var  = data2["r"];
u_var   = data3["u"];
v_var   = data4["v"];

# relative humidity
level=2 # level 2 should correspond to the 700 hPa pressure level. 
lat1 = 1
lat2 = 180

rh_lev       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], dims[3])
u_tot        = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], dims[3])
v_tot        = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], dims[3])
VWS_tot      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], dims[3])

#lon2         = Array{Union{Missing, Float64}, 1}(undef, dims[1])

println("size of rh_lev is: ",size(rh_lev))

println("p level 1 is: ",plev[1])
println("p level 2 is: ",plev[2])

for i in 1:dims[3]
  #rsst_test      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 1)
  rh_lev[:,lat1:lat2,i]   = rh_var[:,lat1:lat2,level,i]
  #rh_lev[:,:,i] = rh_test # ?? should be the average of 6 fields...mpi_b[:,:,]
  u_tot[:,lat1:lat2,i]  = u_var[:,lat1:lat2,2,i] - u_var[:,lat1:lat2,1,i]
  v_tot[:,lat1:lat2,i]  = v_var[:,lat1:lat2,2,i] - v_var[:,lat1:lat2,1,i]
  VWS_tot[:,lat1:lat2,i]  = sqrt.(u_tot[:,lat1:lat2,i].^2 .+ v_tot[:,lat1:lat2,i].^2)
end

#  octNovInd = Matrix{Int32}(undef, 92)
octNovInd = Array{Union{Missing, Float64}, 1}(undef, dims[3])
#pi_var = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], dims[3])
#
## save i and i+1, jump 12, do again
#for i in 10:12:408
#  if i==10 
#    global j=1
#    println("nope!")
#  end
#  octNovInd[j] = i
#  println(j)
#  j = j+1
#  octNovInd[j] = i+1
#  println(j)
#  j = j+1
#end 
#
#print(octNovInd[:])
#pi_var = sst_var[:,:,octNovInd]

print("size of sst Array is: ",size(dim_var))
print("size of PI Array is: ",size(sst_var))
print("size of lon Array is: ",size(lon))
#print(octNovInd[:])

#agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
#bgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])

## loop over longitude and latitude:
#for i in 1:dims[1]
#    for j in 1:dims[2]
#        agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis1[:],pi_var[i,j,:])
#    end
#end
#
#levs = range(-2.0, 2.0, length = 21)
#blah = local_fig(agrid.*20,lon,lat,"PI linear trends October-November (m/s)/decade",levs)
#save("era5_PI_LinTrend_OctThNov_1979th2024_region.png", blah, px_per_unit=6.0)

#save("era5_PI_LinTrend_OctThNov_1979th2024.png", blah, px_per_unit=6.0)
#

#levs = range(-.5, .5, length = 21)

#agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
#for i in 1:dims[1]
#    for j in 1:dims[2]
#        agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis1[:],sst_var[i,j,:])
#    end
#end
#levs = range(-0.5, 0.5, length = 21)
#blah = local_fig(agrid.*20,lon,lat,"sst linear trends October-November (%/decade)",levs)
#save("era5_sst_LinTrend_OctThNov_1979th2024_region.png", blah, px_per_unit=6.0)
##save("era5_sst_LinTrend_OctThNov_1979th2024.png", blah, px_per_unit=6.0)

#levs = range(-0.5, 0.5, length = 21)
#blah = local_fig(agrid.*20,lon,lat,"PI linear trends October-November (m/s)/decade",levs)
#save("era5_PI_LinTrend_OctThNov_1979th2024_region.png", blah, px_per_unit=6.0)


#levs = range(-2.0, 2.0, length = 21)
#blah = local_fig(agrid.*20,lon,lat,"RH linear trends October-November %/decade",levs)
#save("era5_RH850_LinTrend_OctThNov_1979th2024_region.png", blah, px_per_unit=6.0)
##save("era5_RH850_LinTrend_OctThNov_1979th2024.png", blah, px_per_unit=6.0)

# -----------------------------------------
# degrees of freedom 

# average the months for each year to result in arrays with 8 time stamps for 8 years.
# this allows for statistically testing with 7 degrees of freedome instead of 47.

#vws_yrs    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], dims[3]/2)
vws_yrs    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 46)
sst_yrs    = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 46)
pi_yrs     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 46)
rh_yrs     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 46)

jj = 1
#for i in 1:6:48-5
for i in 1:2:dims[3]-1 # assume each year is an independant degree of freedom
    #test  = VWS_tot[:,lat1:lat2,i:i+1];
    test  = VWS_tot[:,:,i:i+1];
    test1 = sst_var[:,:,i:i+1];
    test2 = rh_lev[:,:,i:i+1];
    test3 = pi_var[:,:,i:i+1];
    vws_yrs[:,:,jj]  = mean(test, dims=3)
    sst_yrs[:,:,jj]  = mean(test1, dims=3)
    rh_yrs[:,:,jj]   = mean(test2, dims=3)
    pi_yrs[:,:,jj]   = mean(test3, dims=3)
    global jj += 1
end
# -----------------------------------------

# these arrays store the linear regression values at each grid point
bgrid   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
agrid1  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
agrid2  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
agrid3  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
agrid4  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])

for i in 1:dims[1]
    for j in 1:dims[2]
        agrid1[i,j],bgrid[i,j] = find_best_fit(timeAxis[:],vws_yrs[i,j,:])
        agrid2[i,j],bgrid[i,j] = find_best_fit(timeAxis[:],sst_yrs[i,j,:])
        agrid3[i,j],bgrid[i,j] = find_best_fit(timeAxis[:],rh_yrs[i,j,:])
        #agrid3[i,j],bgrid[i,j] = find_best_fit(timeAxis1[:],rh_lev[i,j,:])
        agrid4[i,j],bgrid[i,j] = find_best_fit(timeAxis[:],pi_yrs[i,j,:])
    end
end

#levs = range(-2.0, 2.0, length = 21)
#blah = local_fig(agrid.*20,lon,lat,"RH linear trends October-November %/decade",levs)
#save("era5_RH850_LinTrend_OctThNov_1979th2024_region.png", blah, px_per_unit=6.0)
##save("era5_RH850_LinTrend_OctThNov_1979th2024.png", blah, px_per_unit=6.0)


##agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
#for i in 1:dims[1]
#    for j in 1:dims[2]
#        #agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis1[:],VWS_tot[i,j,:])
#        agrid1[i,j],bgrid[i,j] = find_best_fit(timeAxis[:],vws_yrs[i,j,:])
#    end
#end

#agrid2 = ShiftedArrays.circshift(agrid, (180,0,0))
#lon2 = collect(0:359.)

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

#println("the arms of the ocean are carrying me ","frac{\sin(x)}{\sqrt{x}}")

    #degree_symbol_string = "25\\u00B0C"
    #degree_symbol_string = "25 \degree"
    #degree_symbol_string = "25" \degree
    #println("New string is: ",degree_symbol_string)



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
    #wilksArray[i] = sortedP[i] - wilks0p05[i]
    wilksArray[i] = sortedP[i] - wilks0p1[i]
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
#println("***********************")

# PI 

pvalgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pvallow   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pvalhigh  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
intVal    = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])

for i in 1:360 # dims[1] longitudes
  for j in 1:180 # dims[2] latitudes
    X2 = Float64.(timeAxis);
    Yind = pi_yrs[i,j,:];
    if any(ismissing, Yind)
      println("missing value found ")
    else
      Ytemp = Float64.(Yind)
      ##
      tab   = Table(X = X2, Y = Ytemp);
      ols_temp = lm(@formula(Y ~ X), tab)
      testVals = coeftable(ols_temp)
      pvalgrid[i,j] = testVals.cols[4][2]
      pvallow[i,j]  = testVals.cols[5][2]
      pvalhigh[i,j] = testVals.cols[6][2]
      intVal[i,j] = sign(pvalhigh[i,j]) + sign(pvallow[i,j])
    end
  end
end

arrN = reshape(pvalgrid, (dims[1]*dims[2]));
sortedP = sort(arrN)

pAxis = collect(1:1:fullength);
wilks0p05 = (0.05/(lengthP)) .* pAxis
wilks0p1 = (0.1/(lengthP)) .* pAxis
p0p05 = zeros(fullength) .+ 0.05;
wilksArray = zeros(fullength);

for i in 1:20000
    #wilksArray[i] = sortedP[i] - wilks0p05[i]
    wilksArray[i] = sortedP[i] - wilks0p1[i]
end

psign = sign.(wilksArray)
diffs = diff(psign)
psign_ind = findall(!iszero, diffs)
println("change of sign should be intersection point of sortedP and wilks0p05: ",psign_ind)
println("PI pvalue at intersection of sorted pvalues and Wilks FDR pvalue is: ",sortedP[psign_ind[1]])

pFDR = sortedP[psign_ind[1]]

# create a bitmatrix to identify statistically significane grid points
woman      = falses(dims[1],dims[2]);
c1 = coalesce.(pvalgrid, false);
c2 = coalesce.(intVal, false);
for i in 1:360
    for j in 1:180
        woman[i,j] = c1[i,j] < pFDR && c2[i,j] != 0.0
    end
end

# grab all the points that are labelled as 'true' or 1.  
points = findall(x -> x == 1, woman); # points has dims of (24065,)
x_coords = [idx.I[1] for idx in points]; # x_coords has dims of (24065,)
y_coords = [idx.I[2] for idx in points];

x_lats_pi = x_coords .- xshift;
y_lons_pi = y_coords .- yshift;

#println("***********************")
#println("***********************")

# RH

pvalgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pvallow   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pvalhigh  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
intVal    = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])

for i in 1:360 # dims[1] longitudes
  for j in 1:180 # dims[2] latitudes
    X2 = Float64.(timeAxis);
    Yind = rh_yrs[i,j,:];
    if any(ismissing, Yind)
      println("missing value found ")
    else
      Ytemp = Float64.(Yind)
      ##
      tab   = Table(X = X2, Y = Ytemp);
      ols_temp = lm(@formula(Y ~ X), tab)
      testVals = coeftable(ols_temp)
      pvalgrid[i,j] = testVals.cols[4][2]
      pvallow[i,j]  = testVals.cols[5][2]
      pvalhigh[i,j] = testVals.cols[6][2]
      intVal[i,j] = sign(pvalhigh[i,j]) + sign(pvallow[i,j])
    end
  end
end

arrN = reshape(pvalgrid, (dims[1]*dims[2]));
sortedP = sort(arrN)

pAxis = collect(1:1:fullength);
wilks0p05 = (0.05/(lengthP)) .* pAxis
wilks0p1 = (0.1/(lengthP)) .* pAxis
p0p05 = zeros(fullength) .+ 0.05;
wilksArray = zeros(fullength);

for i in 1:20000
    #wilksArray[i] = sortedP[i] - wilks0p05[i]
    wilksArray[i] = sortedP[i] - wilks0p1[i]
end

psign = sign.(wilksArray)
diffs = diff(psign)
psign_ind = findall(!iszero, diffs)
println("change of sign should be intersection point of sortedP and wilks0p05: ",psign_ind)
println("RH pvalue at intersection of sorted pvalues and Wilks FDR pvalue is: ",sortedP[psign_ind[1]])

pFDR = sortedP[psign_ind[1]]

# create a bitmatrix to identify statistically significane grid points
woman      = falses(dims[1],dims[2]);
c1 = coalesce.(pvalgrid, false);
c2 = coalesce.(intVal, false);
for i in 1:360
    for j in 1:180
        woman[i,j] = c1[i,j] < pFDR && c2[i,j] != 0.0
    end
end

# grab all the points that are labelled as 'true' or 1.  
points = findall(x -> x == 1, woman); # points has dims of (24065,)
x_coords = [idx.I[1] for idx in points]; # x_coords has dims of (24065,)
y_coords = [idx.I[2] for idx in points];

x_lats_rh = x_coords .- xshift;
y_lons_rh = y_coords .- yshift;

#println("***********************")
#println("***********************")

# VWS 

pvalgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pvallow   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
pvalhigh  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
intVal    = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])

for i in 1:360 # dims[1] longitudes
  for j in 1:180 # dims[2] latitudes
    X2 = Float64.(timeAxis);
    Yind = vws_yrs[i,j,:];
    if any(ismissing, Yind)
      println("missing value found ")
    else
      Ytemp = Float64.(Yind)
      ##
      tab   = Table(X = X2, Y = Ytemp);
      ols_temp = lm(@formula(Y ~ X), tab)
      testVals = coeftable(ols_temp)
      pvalgrid[i,j] = testVals.cols[4][2]
      pvallow[i,j]  = testVals.cols[5][2]
      pvalhigh[i,j] = testVals.cols[6][2]
      intVal[i,j] = sign(pvalhigh[i,j]) + sign(pvallow[i,j])
    end
  end
end

arrN = reshape(pvalgrid, (dims[1]*dims[2]));
sortedP = sort(arrN)

pAxis = collect(1:1:fullength);
wilks0p05 = (0.05/(lengthP)) .* pAxis
wilks0p1 = (0.1/(lengthP)) .* pAxis
p0p05 = zeros(fullength) .+ 0.05;
wilksArray = zeros(fullength);

for i in 1:20000
    #wilksArray[i] = sortedP[i] - wilks0p05[i]
    wilksArray[i] = sortedP[i] - wilks0p1[i]
end

psign = sign.(wilksArray)
diffs = diff(psign)
psign_ind = findall(!iszero, diffs)
println("change of sign should be intersection point of sortedP and wilks0p05: ",psign_ind)
println("VWS pvalue at intersection of sorted pvalues and Wilks FDR pvalue is: ",sortedP[psign_ind[1]])


# the p-values are not well behaved for the VWS so I had to make a few adjustments: 
#pFDR = sortedP[psign_ind[1]]
pFDR = sortedP[psign_ind[2]] # this is for the case of a 0.1 pvalue
#pFDR = sortedP[psign_ind[10]]  # this is best for the case of a 0.05 pvalue

println("VWS pvalue at intersection of sorted pvalues and Wilks FDR pvalue is: ",pFDR)

fig3 = Figure(;
    size = (400,400),
    )
ax = Axis(fig3[1,1];
    xlabel = "p-value rank i",
    ylabel = "p-value",
    title  = "Significance?",
    #limits=(0,2000,0,0.06),
    limits=(0,400,0,0.002),
    )
lines!(ax,pAxis[:],sortedP[:], color = :black)
lines!(ax,pAxis[:],p0p05[:])
lines!(ax,pAxis[:],wilks0p05[:])
lines!(ax,pAxis[:],wilks0p1[:], color = :red)
fig3
save("sig_Wilks_vws.png", fig3, px_per_unit=6.0)

#----------------------
# create a bitmatrix to identify statistically significane grid points
woman      = falses(dims[1],dims[2]);
c1 = coalesce.(pvalgrid, false);
c2 = coalesce.(intVal, false);
for i in 1:360
    for j in 1:180
        woman[i,j] = c1[i,j] < pFDR && c2[i,j] != 0.0
        #woman[i,j] = c1[i,j] < 0.05 && c2[i,j] != 0.0
    end
end

# grab all the points that are labelled as 'true' or 1.  
points = findall(x -> x == 1, woman); # points has dims of (24065,)
x_coords = [idx.I[1] for idx in points]; # x_coords has dims of (24065,)
y_coords = [idx.I[2] for idx in points];

x_lats_vws = x_coords .- xshift;
y_lons_vws = y_coords .- yshift;

#println("***********************")
#println("***********************")




# end of statistical analysis
#println("***********************")

# print figure panels

#levs = range(-2.0, 2.0, length = 21)
##blah = local_fig(agrid.*20,lon,lat,"VWS linear trends October-November (m/s)/decade",levs)
#blah = local_fig(agrid.*20,lon,lat,"VWS linear trends October-November (m/s)/decade",levs,:vik)
#save("era5_VWS_LinTrend_OctThNov_1979th2024_region.png", blah, px_per_unit=6.0)
##save("era5_VWS_LinTrend_OctThNov_1979th2024.png", blah, px_per_unit=6.0)

levs = range(-5.0, 5.0, length = 21)
blah = fig_plot_stats(agrid3.*20,lon,lat,"700 hPa Relative Humidity linear trends October-November (%)/decade",levs,:BrBg, x_lats_rh, y_lons_rh)
save("era5_RH_LinTrend_OctThNov_1979th2024_region_46dof.png", blah, px_per_unit=6.0)

## sst levels
levs = range(-0.5, 0.5, length = 21)
blah = fig_plot_stats(agrid2.*20,lon,lat,"SST linear trends October-November ( C)/decade",levs,:vik, x_lats_sst, y_lons_sst)
save("era5_SST_LinTrend_OctThNov_1979th2024_region_46dof.png", blah, px_per_unit=6.0)
#
## VWS levels
levs = range(-3.0, 3.0, length = 21)
blah = fig_plot_stats(agrid1.*20,lon,lat,"Vertical Wind Shear linear trends October-November (m/s)/decade",levs,:bam, x_lats_vws, y_lons_vws)
save("era5_VWS_LinTrend_OctThNov_1979th2024_region_46dof.png", blah, px_per_unit=6.0)

## PI levels
levs = range(-5.0, 5.0, length = 21)
blah = fig_plot_stats(agrid4.*20,lon,lat,"Potential Intensity linear trends October-November (m/s)/decade",levs,:vik, x_lats_pi, y_lons_pi)
save("era5_PI_LinTrend_OctThNov_1979th2024_region_46dof.png", blah, px_per_unit=6.0)

#blah = fig_plot_stats(agrid1.*20,lon,lat,"VWS linear trends October-November (m/s)/decade",levs,:bam, x_lats, y_lons)
#save("era5_VWS_LinTrend_OctThNov_1979th2024_region_46dof.png", blah, px_per_unit=6.0)






