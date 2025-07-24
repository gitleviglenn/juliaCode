#---------------------------------------------------------------------
# RelativeSST_PI_RH_VWS_Corr.jl
#
# levi silvers                                               july 2025
#---------------------------------------------------------------------

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics
using DataFrames
using CSV

include("ensoFuncs.jl")

path="/Users/C823281551/data/ERA5/"
path2 = "/Users/C823281551/data/fromExcel/"

filein  = path*"era5_sst_1990th2023_360x180.nc"
file2   = path*"MPI_ERA5_full_output.nc"
file3   = path2*"NA_RI_storms_1990th2023.csv"
tag = "ERA5"
data   = NCDataset(filein);
data2  = NCDataset(file2);
data3  = CSV.read(file3, DataFrame)

lat = data["lat"]
lon = data["lon"]
tme = data["valid_time"]
timeAxis1 = collect(1.083333:1/12:35);
timeAxis2 = collect(1.083333:1/6:35);
timeAxis3 = collect(1.:1:204);
timeENSO  = collect(1.:1:48);

sst_var    = data["sst"];
vmax       = data2["vmax"];

ri_events  = data3.NumStorms

function create_indices(years)
  ensoInd = Matrix{Int32}(undef, 34, 6)
  for i in 1:34
    si = years[i]
    ensoInd[i,:] = si:si+5
  end
  return ensoInd
end

# create an array for the initial index of June for each year
yearsArray = collect(6:12:408)

# create an array with 6 columns (months) for each year
yearsInd = create_indices(yearsArray[:])

# define our region of interest to be +/- 40?
lat1 = 50
lat2 = 130

latS = 60
latN = 120
endIndex = 204

dims = size(sst_var)
sst_trm      = Array{Union{Missing, Float64}, 1}(undef, endIndex)
#sst_trm     = Array{Union{Missing, Float64}, 3}(undef, dims[1], 80, 6)
sst_rel      = Array{Union{Missing, Float64}, 3}(undef, dims[1], 81, endIndex)
mpi          = Array{Union{Missing, Float64}, 3}(undef, dims[1], 81, endIndex)

print("size of yearsArray is: ",size(yearsArray))
print("size of yearsInd is: ",size(yearsInd))
print("size of full sst_var array is: ",size(sst_var))

print("size of ri_events: ", ri_events)

println("yearsInd is: ",yearsInd)
# compute tropical mean
for i in 1:endIndex
    #sst_trm_nino[i] = mean(skipmissing(sst_nino[:,lat1:lat2,i]))
    #sst_trm_nina[i] = mean(skipmissing(sst_nina[:,lat1:lat2,i]))
    #sst_trm_nina[i] = mean(skipmissing(sst_var[:,latS:latN,ninaInd[i]]))
    sst_trm[i] = mean(skipmissing(sst_var[:,latS:latN,yearsInd[i]]))
end

for i in 1:endIndex
  #sst_rel[:,lat1:lat2,i]  = sst_var[:,lat1:lat2,yearsInd[i]] .- sst_trm[i]
  #mpi[:,lat1:lat2,i]      = vmax[:,lat1:lat2,yearsInd[i]]
  sst_rel[:,:,i]  = sst_var[:,lat1:lat2,yearsInd[i]] .- sst_trm[i]
  mpi[:,:,i]      = vmax[:,lat1:lat2,yearsInd[i]]
end

dims = size(sst_var)

#mpi_nino          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
mpi_yr     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 34)
mpi_1yr    = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])


# another way...
sst_trm_b      = Array{Union{Missing, Float64}, 1}(undef, endIndex)
sst_rel_b      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], endIndex)
mpi_b          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 6)

for i in 1:34
  mpi_test       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 1)
  ii = i
  for j in 1:6 # compute tropical mean for each month of interest
    sst_trm_b[j] = mean(skipmissing(sst_var[:,latS:latN,yearsInd[i,j]]))
  end
  for j in 1:6
    sst_rel_b[:,lat1:lat2,j]  = sst_var[:,lat1:lat2,yearsInd[i,j]] .- sst_trm_b[j]
    mpi_b[:,lat1:lat2,j]      = vmax[:,lat1:lat2,yearsInd[i,j]]
    #println("j value is: ",j,", and yearsInd is: ",yearsInd[i,j])
  end
  #println("ii value is: ",ii)
  mpi_test = mean(mpi_b, dims=3)
  mpi_yr[:,:,i] = mpi_test #should be the average of 6 fields...mpi_b[:,:,]
end

levs = range(60., 90., length = 21)
# colormap = :vik, seems to work well for MPI
#mpiDiff = fig_anom_plot(mpi_b[:,:,1],lon,lat,"MPI (m/s), index 1",levs)
mpiDiff = fig_1_plot(mpi_b[:,:,1],lon,lat,"MPI (m/s), index 1",levs)
save("era5_mpi_1.png", mpiDiff, px_per_unit=6.0)

mpiDiff = fig_1_plot(mpi_b[:,:,2],lon,lat,"MPI (m/s), index 1",levs)
save("era5_mpi_2.png", mpiDiff, px_per_unit=6.0)
#
mpiDiff = fig_1_plot(mpi_b[:,:,3],lon,lat,"MPI (m/s), index 1",levs)
save("era5_mpi_3.png", mpiDiff, px_per_unit=6.0)
#
mpiDiff = fig_1_plot(mpi_b[:,:,4],lon,lat,"MPI (m/s), index 1",levs)
save("era5_mpi_4.png", mpiDiff, px_per_unit=6.0)
#
mpiDiff = fig_1_plot(mpi_b[:,:,5],lon,lat,"MPI (m/s), index 1",levs)
save("era5_mpi_5.png", mpiDiff, px_per_unit=6.0)
#
mpiDiff = fig_1_plot(mpi_b[:,:,6],lon,lat,"MPI (m/s), index 1",levs)
save("era5_mpi_6.png", mpiDiff, px_per_unit=6.0)

mpiDiff = fig_1_plot(mpi_yr[:,:,1],lon,lat,"MPI (m/s), index 1",levs)
#mpiDiff = fig_1_plot(mpi_test[:,:,1],lon,lat,"MPI (m/s), index 1",levs)
save("era5_mpi_yr_1.png", mpiDiff, px_per_unit=6.0)

mpiDiff = fig_1_plot(mpi_yr[:,:,2],lon,lat,"MPI (m/s), index 1",levs)
#mpiDiff = fig_1_plot(mpi_test[:,:,1],lon,lat,"MPI (m/s), index 1",levs)
save("era5_mpi_yr_2.png", mpiDiff, px_per_unit=6.0)

mpiDiff = fig_1_plot(mpi_yr[:,:,12],lon,lat,"MPI (m/s), index 1",levs)
save("era5_mpi_yr_12.png", mpiDiff, px_per_unit=6.0)










