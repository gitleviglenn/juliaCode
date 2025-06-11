#-------------------------------------------------------------------------------
# RelativeSST_PI_ENSO.jl
#
# Compute the relative SST change by using a linear regression calculation
# data: ERA5
#
# Default behavior is not to compute the linear regression sst change of each 
# grid point to change this set regThres > 2
#
# levi silvers                                                    June 2025
#-------------------------------------------------------------------------------

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics
#using GLM
using DataFrames

include("ensoFuncs.jl")

path="/Users/C823281551/data/ERA5/"

filein  = path*"era5_sst_1990th2023_360x180.nc"
file2   = path*"era5_sst_1990th2023_360x180_summer.nc"
file3   = path*"MPI_ERA5_full_output.nc"
tag = "ERA5"
data   = NCDataset(filein);
data2  = NCDataset(file2);
data3  = NCDataset(file3);

lat = data["lat"]
lon = data["lon"]
tme = data["valid_time"]
timeAxis1 = collect(1.083333:1/12:35);
timeAxis2 = collect(1.083333:1/6:35);
timeAxis3 = collect(1.:1:204);
timeENSO  = collect(1.:1:48);

sst_var    = data["sst"];
sst_summer = data2["sst"];
vmax       = data3["vmax"];

# Create index lists for ENSO years, and grab the corresponding data from sst_var[]

## Northern Hemisphere
ninoyears = [18 54 90 150 174 234 306 402]
ninayears = [102 114 210 246 318 366 378 390]
# Southern Hemisphere
# ninoyears = [23 35 59 96 155 239 311 347]
# ninayears = [107 119 215 251 263 335 371 383]

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

dims = size(sst_var)
dimh = Int.(dims[3]/2)
numfields = 48
ensoEnd = 48

sst_nino          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
sst_nina          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
mpi_nino          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
mpi_nina          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
sst_trm_nino      = Array{Union{Missing, Float64}, 1}(undef, ensoEnd)
sst_trm_nina      = Array{Union{Missing, Float64}, 1}(undef, ensoEnd)

# define our region of interest to be +/- 40?
lat1 = 50
lat2 = 130

latS = 60
latN = 120

for i in 1:ensoEnd
    #sst_trm_nino[i] = mean(skipmissing(sst_nino[:,lat1:lat2,i]))
    #sst_trm_nina[i] = mean(skipmissing(sst_nina[:,lat1:lat2,i]))
    sst_trm_nina[i] = mean(skipmissing(sst_var[:,latS:latN,ninaInd[i]]))
    sst_trm_nino[i] = mean(skipmissing(sst_var[:,latS:latN,ninoInd[i]]))
end

for i in 1:48
  sst_nino[:,lat1:lat2,i] = sst_var[:,lat1:lat2,ninoInd[i]] .- sst_trm_nino[i]
  sst_nina[:,lat1:lat2,i] = sst_var[:,lat1:lat2,ninaInd[i]] .- sst_trm_nina[i]
  mpi_nino[:,lat1:lat2,i] = vmax[:,lat1:lat2,ninoInd[i]]
  mpi_nina[:,lat1:lat2,i] = vmax[:,lat1:lat2,ninaInd[i]]
end

# Compute time averages

#dimsNino = size(sst_nino)
#ensoEnd = dimsNino[3]

#sst_tr_mean_full     = Array{Union{Missing, Float64}, 1}(undef, 408)
#sst_tr_mean          = Array{Union{Missing, Float64}, 1}(undef, 204)
#sst_mn = mean(sst_var, dims=3)
#sst_tr_mn = mean.(skipmissing(sst_var[:,lat1:lat2,:]));

# i don't think we need these time averages any more...
#for i in 1:408
#    sst_tr_mean_full[i] = mean(skipmissing(sst_var[:,lat1:lat2,i]))
#end
#for i in 1:204
#    sst_tr_mean[i] = mean(skipmissing(sst_summer[:,lat1:lat2,i]))
#end


#nino_mn = mean(sst_nino, dims=3)
#nina_mn = mean(sst_nina, dims=3)
#mpi_nino_mn = mean(mpi_nino, dims=3)
#mpi_nina_mn = mean(mpi_nina, dims=3)

print("size of MPI Nina mean is: ",size(mpi_nina))

mpi_comp    = mpi_nino .- mpi_nina
mpi_comp_mn = mean(mpi_comp, dims=3)

sst_comp    = sst_nino .- sst_nina
sst_comp_mn = mean(sst_comp, dims=3)

# begin work on the relative SST composite figure: 
#--------------------------------------------------------------------------------------------------------------

# if linear regression is desired as the means of computing relative sst change at a 
# grid point then set regThres > 1

regThres = 0

if regThres > 1
    # use linear regression to model the trend at each grid point.
    agrid1  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
    bgrid1  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
    anino   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
    bnino   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
    anina   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
    bnina   = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
    # loop over longitude and latitude:
    for i in 1:dims[1]
        for j in 1:dims[2]
            agrid1[i,j],bgrid1[i,j] = find_best_fit(timeAxis3[:],sst_summer[i,j,:])
            ##agrid1[i,j],bgrid1[i,j] = find_best_fit(timeENSO[:],sst_summer[i,j,:])
            ##agrid1[i,j],bgrid1[i,j] = find_best_fit(timeAxis2[:],sst_summer[i,j,:].-297.03)
            #anino[i,j],bnino[i,j] = find_best_fit(timeENSO[:],sst_nino[i,j,:])
            #anina[i,j],bnina[i,j] = find_best_fit(timeENSO[:],sst_nina[i,j,:])
        end
    end
    # timeAxis2 consists of monthly values over 204 timesteps taken from 34 years.  So our trend values should be
    # Degrees/34 years right?  Then to convert to degrees per century 'agrid1' could be multiplied by # 2.94
    
    
    
    
    
    #levs = range(-2., 2., length = 21)
    #blah = fig_1_anom(agrid1.*600,lon,lat,levs,"linear trends C per Century")
    ## is this the equivalent of Figure 3a in Vecchi and Soden, 2007 but for ERA5 over a different time period?
    #save("era5_test.png", blah, px_per_unit=6.0)
    
    # subtract off tropical mean.
    relSST = agrid1*600 .- 1.55; # 1.55 C/century, estimated from the slope of the linear trend of ERA5 tropical mean.

end

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

levs = range(-10., 10., length = 21)
#mpiDiff = fig_1_anom(mpi_comp[:,:,1],lon,lat,levs,"Composite MPI (m/s), El Nino - La Nina")

# colormap = :vik, seems to work well for MPI
mpiDiff = fig_anom_plot(mpi_comp_mn[:,:,1],lon,lat,"Composite MPI (m/s), El Nino - La Nina",levs)
#mpiDiff = fig_anom_plot(mpi_comp[:,:,2],lon,lat,"Composite MPI (m/s), El Nino - La Nina",levs)
save("era5_mpiDiff_NH.png", mpiDiff, px_per_unit=6.0)

levs = range(-2., 2., length = 21)
sstDiff = fig_anom_plot(sst_comp_mn[:,:,1],lon,lat,"Composite MPI (m/s), El Nino - La Nina",levs)
save("era5_sstDiff_NH.png", sstDiff, px_per_unit=6.0)







