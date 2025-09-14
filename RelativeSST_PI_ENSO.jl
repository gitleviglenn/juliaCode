#-------------------------------------------------------------------------------
# RelativeSST_PI_ENSO.jl
#
# Compute the relative SST change by using a linear regression calculation
# data: ERA5
#
# levi silvers                                                    June 2025
#-------------------------------------------------------------------------------

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics
using HypothesisTests
using DataFrames
using TypedTables

include("ensoFuncs.jl")

path="/Users/C823281551/data/ERA5/"

filein  = path*"era5_sst_1990th2023_360x180.nc"
#file2   = path*"era5_sst_1990th2023_360x180_summer.nc"
file3   = path*"MPI_ERA5_full_output.nc"
tag = "ERA5"
data   = NCDataset(filein);
#data2  = NCDataset(file2);
data3  = NCDataset(file3);

lat = data["lat"]
lon = data["lon"]
tme = data["valid_time"]
#timeAxis1 = collect(1.083333:1/12:35);
timeAxis2 = collect(1.083333:1/6:35);
#timeAxis3 = collect(1.:1:204);
#timeAxis4 = collect(1.:1:408);
timeENSO  = collect(1.:1:48);

sst_var    = data["sst"];
#sst_summer = data2["sst"];
vmax       = data3["vmax"];

# Create index lists for ENSO years, and grab the corresponding data from sst_var[]

## Northern Hemisphere
ninoyears = [18 54 90 150 174 234 306 402]
ninayears = [102 114 210 246 318 366 378 390]
## Southern Hemisphere
#ninoyears = [23 35 59 96 155 239 311 347]
#ninayears = [107 119 215 251 263 335 371 383]

#----------
# Parameters for SST: 
#figname = "era5_sstDiff_SH_statWilks_8dof.png"
#maintit="Relative SST Composite" 
#levs = range(-2., 2., length = 21)
#----------
# Parameters for MPI: 
figname = "era5_MPI_NH_statWilks_8dof.png"
maintit="Potential Intensity Composite"
levs = range(-10., 10., length = 21)

function create_indices(years)
  ensoInd = Matrix{Int32}(undef, 8, 6)
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

println("-----------------------------------------------------------------------------------")
println("Nino indices are: ",ninoInd[1,:])
println("Nino indices are: ",ninoInd[2,:])
println("Nino indices are: ",ninoInd[3,:])
println("Nino indices are: ",ninoInd[4,:])
println("Nino indices are: ",ninoInd[5,:])
println("Nino indices are: ",ninoInd[6,:])
println("Nino indices are: ",ninoInd[7,:])
println("Nino indices are: ",ninoInd[8,:])
println("-----------------------------------------------------------------------------------")
println("Nina indices are: ",ninaInd[1,:])
println("Nina indices are: ",ninaInd[2,:])
println("Nina indices are: ",ninaInd[3,:])
println("Nina indices are: ",ninaInd[4,:])
println("Nina indices are: ",ninaInd[5,:])
println("Nina indices are: ",ninaInd[6,:])
println("Nina indices are: ",ninaInd[7,:])
println("Nina indices are: ",ninaInd[8,:])
println("-----------------------------------------------------------------------------------")

dims = size(sst_var)
println("size of sst_var is: ",size(sst_var))
dimh = Int.(dims[3]/2)
numfields = 48
ensoEnd = 48

sst_nino          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
sst_nina          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
mpi_nino          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
mpi_nina          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
sst_trm_nino      = Array{Union{Missing, Float64}, 1}(undef, ensoEnd)
sst_trm_nina      = Array{Union{Missing, Float64}, 1}(undef, ensoEnd)
sst_enso_yrs_high = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 8)
sst_enso_yrs_low  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 8)
mpi_enso_yrs_high = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 8)
mpi_enso_yrs_low  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 8)

# define our region of interest to be +/- 40?
lat1 = 50
lat2 = 130

latS = 60
latN = 120

for i in 1:ensoEnd
  sst_trm_nina[i] = mean(skipmissing(sst_var[:,latS:latN,ninaInd[i]]))
  sst_trm_nino[i] = mean(skipmissing(sst_var[:,latS:latN,ninoInd[i]]))
end

for i in 1:48
  sst_nino[:,lat1:lat2,i] = sst_var[:,lat1:lat2,ninoInd[i]] .- sst_trm_nino[i]
  sst_nina[:,lat1:lat2,i] = sst_var[:,lat1:lat2,ninaInd[i]] .- sst_trm_nina[i]
  mpi_nino[:,lat1:lat2,i] = vmax[:,lat1:lat2,ninoInd[i]]
  mpi_nina[:,lat1:lat2,i] = vmax[:,lat1:lat2,ninaInd[i]]
end

# average the months for each year to result in arrays with 8 time stamps for 8 years.
# this allows for statistically testing with 7 degrees of freedome instead of 47.
jj = 1
for i in 1:6:48-5
    #println("i index is: ",i)
    #println("j index is: ",j)
    test  = sst_nina[:,lat1:lat2,i:i+5];
    test2 = sst_nino[:,lat1:lat2,i:i+5];
    sst_enso_yrs_low[:,lat1:lat2,jj]  = mean(test, dims=3)
    sst_enso_yrs_high[:,lat1:lat2,jj] = mean(test2, dims=3)
    test3 = mpi_nina[:,lat1:lat2,i:i+5];
    test4 = mpi_nino[:,lat1:lat2,i:i+5];
    mpi_enso_yrs_low[:,lat1:lat2,jj]  = mean(test3, dims=3)
    mpi_enso_yrs_high[:,lat1:lat2,jj] = mean(test4, dims=3)
    global jj += 1
end

# Compute time averages

#dimsNino = size(sst_nino)
#ensoEnd = dimsNino[3]

#sst_tr_mean_full     = Array{Union{Missing, Float64}, 1}(undef, 408)
#sst_tr_mean          = Array{Union{Missing, Float64}, 1}(undef, 204)
#sst_mn = mean(sst_var, dims=3)
#sst_tr_mn = mean.(skipmissing(sst_var[:,lat1:lat2,:]));

#nino_mn = mean(sst_nino, dims=3)
#nina_mn = mean(sst_nina, dims=3)
#mpi_nino_mn = mean(mpi_nino, dims=3)
#mpi_nina_mn = mean(mpi_nina, dims=3)

println("size of full vmax array is: ",size(vmax))
println("size of MPI Nina mean is: ",size(mpi_nina))

mpi_comp    = mpi_nino .- mpi_nina
mpi_comp_mn = mean(mpi_comp, dims=3)

sst_comp    = sst_nino .- sst_nina
sst_comp_mn = mean(sst_comp, dims=3)

sst_comp_8dof = sst_enso_yrs_high .- sst_enso_yrs_low
mpi_comp_8dof = mpi_enso_yrs_high .- mpi_enso_yrs_low

println("size of sst_comp is: ",size(sst_comp))
println("size of sst_comp_mn is: ",size(sst_comp_mn))

# begin work on the relative SST composite figure: 
#--------------------------------------------------------------------------------------------------------------

#  # use linear regression to model the trend at each grid point.
#  agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
#  bgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
#endTime = 60
#    for i in 1:dims[1]
#        for j in 1:dims[2]
#   ##         agrid[i,j],bgrid[i,j] = find_best_fit(timeENSO[:],sst_comp[i,j,:])
#   #         agrid[i,j],bgrid[i,j] = find_best_fit(timeAxis4[1:endTime],sst_var[i,j,1:endTime])
#        end
#    end
#    # timeAxis2 consists of monthly values over 204 timesteps taken from 34 years.  So our trend values should be
#    # Degrees/34 years right?  Then to convert to degrees per century 'agrid1' could be multiplied by # 2.94
#    #println("agrid: ",agrid[50:100,50:100])
#    levs = range(-0.1, 0.1, length = 21)
#    blah = fig_anom_plot(agrid,lon,lat,"linear trends C per Century",levs)
#    ## is this the equivalent of Figure 3a in Vecchi and Soden, 2007 but for ERA5 over a different time period?
#    save("era5_sstlintrend.png", blah, px_per_unit=6.0)
    
    # subtract off tropical mean.
#    relSST = agrid1*600 .- 1.55; # 1.55 C/century, estimated from the slope of the linear trend of ERA5 tropical mean.

#end

#X2 = Float64.(timeENSO);
#X2 = Float64.(timeAxis4[1:endTime]);

d1 = 360
d2 = 81
pvalgrid  = Array{Union{Missing, Float64}, 2}(undef, d1, d2)
pvallow   = Array{Union{Missing, Float64}, 2}(undef, d1, d2)
pvalhigh  = Array{Union{Missing, Float64}, 2}(undef, d1, d2)
intVal    = Array{Union{Missing, Float64}, 2}(undef, d1, d2)

for i in 1:360 # dims[1] longitudes
    #for j in 1:180 # dims[2] latitudes
    for j in 50:130 # dims[2] latitudes
        #Yind = sst_comp[i,j,:];
        #Yind = sst_comp_8dof[i,j,:];
        #Yind = mpi_comp[i,j,:];
        Yind = mpi_comp_8dof[i,j,:];
        if any(ismissing, Yind)
            #println("missing value found ")
        else
        # need to somehow check if Yind contains missing data, and if so, then don't proceed to the lm step.  
          Ytemp = Float64.(Yind)
          ols_temp = OneSampleTTest(Ytemp)
          pvalgrid[i,j-49] = pvalue(ols_temp)
          rang = confint(ols_temp, level = 0.95, tail = :both);
          pvallow[i,j-49]  = rang[1]
          pvalhigh[i,j-49] = rang[2]
          intVal[i,j-49] = sign(pvalhigh[i,j-49]) + sign(pvallow[i,j-49])
        end
    end
end

arrN = reshape(pvalgrid, (d1*d2));
sortedP = sort(arrN)

lengthP = 81*360;
#arrN = reshape(pvalgrid, (lengthP));
#testSort = sort(arrN)

pAxis = collect(1:1:lengthP);
wilks0p05 = (0.05/(lengthP)) .* pAxis
wilks0p1 = (0.1/(lengthP)) .* pAxis
p0p05 = zeros(lengthP) .+ 0.05;
wilksArray = zeros(lengthP);
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

#blab = where(testSort = wilks0p05)
# blah = maximum (testSort <= wilks0p05)
for i in 1:20000
    #wilksArray[i] = testSort[i] <= wilks0p05[i]
    wilksArray[i] = sortedP[i] - wilks0p05[i]
end
#println(wilks0p05[1:500])
#println("*******************")
#println(testSort[1:500])
#println("wilks threshold: ",maximum(wilksArray))

#lines!(ax,pAxis[:],wilks0p1[:])
fig3
save("sig_Wilks.png", fig3, px_per_unit=6.0)

fig3 = Figure(;
    size = (400,400),
    )
ax = Axis(fig3[1,1];
    xlabel = "p-value rank i",
    ylabel = "p-value",
    title  = "Significance?",
    limits=(0,20000,-0.1,0.1),
    )
#lines!(C,ts_roni_sm, linestyle = :solid)
lines!(ax,pAxis[:],wilksArray[:], color = :black)
fig3
save("sig_WilksArray.png", fig3, px_per_unit=6.0)


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
    #for j in 50:130
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
x_lats = x_coords .- 181;
#y_lons = y_coords .- 90;
y_lons = y_coords .- 42;

#levs = range(-2., 2., length = 21)
#tropTrend = fig_1_anom(relSST,lon,lat,levs,"linear trends C per century, minus trop mean")
## is this the equivalent of Figure 3b in Vecchi and Soden, 2007 but for ERA5 over a different time period?
#save("era5_TropTrend_minusMean.png", tropTrend, px_per_unit=6.0)

# what do I want to show in my paper?   Simply the composite values of SST and MPI?   Or perhaps the relative MPI?

#levs = range(40., 100., length = 61)
##mpiNina = contourf(mpi_nina_mn[:,:,1])
#mpiNina = fig_1_plot(mpi_nina_mn[:,:,1],lon,lat,levs,"MPI during La Nina")
## is this the equivalent of Figure 3b in Vecchi and Soden, 2007 but for ERA5 over a different time period?
#save("era5_mpi_nina.png", mpiNina, px_per_unit=6.0)

#levs = range(-10., 10., length = 21)
##mpiDiff = fig_1_anom(mpi_comp[:,:,1],lon,lat,levs,"Composite MPI (m/s), El Nino - La Nina")
#
## colormap = :vik, seems to work well for MPI
#mpiDiff = fig_anom_plot(mpi_comp_mn[:,:,1],lon,lat,"Composite MPI (m/s), El Nino - La Nina",levs)
##mpiDiff = fig_anom_plot(mpi_comp[:,:,2],lon,lat,"Composite MPI (m/s), El Nino - La Nina",levs)
#save("era5_mpiDiff_NHb.png", mpiDiff, px_per_unit=6.0)

#levs = range(-2., 2., length = 21)
#sstDiff = fig_anom_plot(sst_comp_mn[:,:,1],lon,lat,"Composite SST, El Nino - La Nina",levs)
#save("era5_sstDiff_NHb.png", sstDiff, px_per_unit=6.0)



# range for SST: 
#levs = range(-2., 2., length = 21)
# range for MPI: 
#levs = range(-10., 10., length = 21)
    f3 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:white,                                                                     
        size=(900,400),                                                                             
        )   
    ax = GeoAxis(f3[1,1];
        #aspect = 3,
        #dest="+proj=latlon",                                                                        
        xticks = -180:30:180,                                                                       
        yticks = -90:30:90,                                                                         
        xlabel="longitude",
        ylabel="latitude",
        limits=(-180,180,-40,40),
        #title="Relative SST Composite", 
        #title="Potential Intensity Composite",
        title=maintit, 
        xticklabelsize = 22, # 14,16 are pretty reasonable sizes                                    
        yticklabelsize = 22, # 22 used for 8 panel figure that needs larger font                    
        )                                                                                           
        #bb = contourf!(ax, lon, lat, sst_comp_mn[:,:,1],                                                            
        bb = contourf!(ax, lon, lat, mpi_comp_mn[:,:,1],                                                            
             levels = levs,
             #colormap = :batlow,
             #colormap = :bam, # default for shear plot (greens and pinks)
             #colormap = :seismic, # colors are a bit harsh                                         
             colormap = :vik, # default for redish bluish for relative SST                          
             #colormap = :BrBg, # better for RH  browns and greens                                  
             #colormap = :roma,
             extendlow = :auto, extendhigh = :auto                                                  
        )       
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)                           
        Colorbar(f3[1,2], bb)                                                                       
scatter!(ax, x_lats, y_lons, marker = :circle, markersize=2, color = :black)
f3

save(figname, f3, px_per_unit=6.0)
#save("era5_sstDiff_SH_statWilks_8dof.png", f3, px_per_unit=6.0)

#levs = range(-0.1, 0.1, length = 21)
#    f3 = Figure(;
#        figure_padding=(5,5,10,10),
#        backgroundcolor=:white,                                                                     
#        size=(900,400),                                                                             
#        )   
#    ax = GeoAxis(f3[1,1];
#        aspect = 3,
#        dest="+proj=latlon", # what projections corresponds to 1:1?
#        xticks = -180:20:180,                                                                       
#        yticks = -40:20:40,                                                                         
#        xlabel="longitude",
#        ylabel="latitude",
#        limits=(-180,180,-40,40),
#        title="Relative SST Composite", 
#        #title="Potential Intensity Composite", 
#        xticklabelsize = 14, # 14,16 are pretty reasonable sizes                                    
#        yticklabelsize = 14, # 22 used for 8 panel figure that needs larger font                    
#        )                                                                                           
#        #bb = contourf!(ax, lon, lat, sst_comp_mn[:,:,1],                                                            
#        bb = contourf!(ax, lon, lat, agrid[:,:],                                                            
#             levels = levs,
#             #colormap = :batlow,
#             #colormap = :bam, # default for shear plot (greens and pinks)
#             #colormap = :seismic, # colors are a bit harsh                                         
#             colormap = :vik, # default for redish bluish for relative SST                          
#             #colormap = :BrBg, # better for RH  browns and greens                                  
#             #colormap = :roma,
#             extendlow = :auto, extendhigh = :auto                                                  
#        )       
#        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)                           
#        Colorbar(f3[1,2], bb)                                                                       
#scatter!(ax, x_lats, y_lons, marker = :utriangle, markersize=4, color = :black)
#f3
#
#save("era5_sstLinTrend_statWilks.png", f3, px_per_unit=6.0)


