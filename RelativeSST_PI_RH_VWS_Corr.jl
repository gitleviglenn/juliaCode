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
file3   = path*"era5_rh_1990to2023_360x180.nc"
file4  = path*"era5_uv_1990th2023_360x180.nc"
#file3   = path*"era5_rh_1990th2023_360x180.nc"
file5   = path2*"NA_RI_storms_1990th2023.csv"
tag = "ERA5"
data   = NCDataset(filein);
data2  = NCDataset(file2);
data3  = NCDataset(file3);
data4  = NCDataset(file4);
data5  = CSV.read(file5, DataFrame)

lat = data["lat"]
lon = data["lon"]
lev = data3["pressure_level"]
tme = data["valid_time"]
timeAxis1 = collect(1.083333:1/12:35);
timeAxis2 = collect(1.083333:1/6:35);
timeAxis3 = collect(1.:1:204);
timeENSO  = collect(1.:1:48);

# read variables into arrayes/dataframes:
sst_var    = data["sst"];
vmax       = data2["vmax"];
rh_var     = data3["r"];
u_var      = data4["u"];
v_var      = data4["v"];
ri_events  = data5.NumStorms

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
rh           = Array{Union{Missing, Float64}, 3}(undef, dims[1], 81, endIndex)

print("size of yearsArray is: ",size(yearsArray))
print("size of yearsInd is: ",size(yearsInd))
print("size of full sst_var array is: ",size(sst_var))

print("size of ri_events: ", ri_events)

println("yearsInd is: ",yearsInd)
## compute tropical mean
#for i in 1:endIndex
#    sst_trm[i] = mean(skipmissing(sst_var[:,latS:latN,yearsInd[i]]))
#end
#
#for i in 1:endIndex
#  #sst_rel[:,lat1:lat2,i]  = sst_var[:,lat1:lat2,yearsInd[i]] .- sst_trm[i]
#  #mpi[:,lat1:lat2,i]      = vmax[:,lat1:lat2,yearsInd[i]]
#  sst_rel[:,:,i]  = sst_var[:,lat1:lat2,yearsInd[i]] .- sst_trm[i]
#  mpi[:,:,i]      = vmax[:,lat1:lat2,yearsInd[i]]
#  rh[:,:,i]       = rh_var[:,lat1:lat2,level,yearsInd[i]]
#end

dims = size(sst_var)

mpi_yr      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 34)
rsst_yr     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 34)
rh_yr       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 34)
vws_yr      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 34)

# another way...
sst_trm_b      = Array{Union{Missing, Float64}, 1}(undef, endIndex)
sst_rel_b      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 6)
mpi_b          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 6)
rh_b           = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 6)
u_tot_b        = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 6)
v_tot_b        = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 6)
VWS_tot_b      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 6)

level=2 # level 2 should correspond to the 700 hPa pressure level. 

println("------------------------------------")
print("size of vmax is: ",size(vmax))
print("size of rh_var is: ",size(rh_var))
print("size of rh is: ",size(rh))
println("------------------------------------")

#for i in 1:408
#  # calculate the shear for all times
#  u_tot[:,lat1:lat2,i]   = u_var[:,lat1:lat2,2,i] - u_var[:,lat1:lat2,1,i]
#  v_tot[:,lat1:lat2,i]   = v_var[:,lat1:lat2,2,i] - v_var[:,lat1:lat2,1,i]
#  VWS_tot[:,lat1:lat2,i] = sqrt.(u_tot[:,lat1:lat2,i].^2 .+ v_tot[:,lat1:lat2,i].^2)
#end

for i in 1:34
  mpi_test       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 1)
  rsst_test      = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 1)
  rh_test        = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 1)
  vws_test       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], 1)
  ii = i
  for j in 1:6 # compute tropical mean for each month of interest
    sst_trm_b[j] = mean(skipmissing(sst_var[:,latS:latN,yearsInd[i,j]]))
  end
  for j in 1:6
    sst_rel_b[:,lat1:lat2,j]  = sst_var[:,lat1:lat2,yearsInd[i,j]] .- sst_trm_b[j]
    mpi_b[:,lat1:lat2,j]      = vmax[:,lat1:lat2,yearsInd[i,j]]
    rh_b[:,lat1:lat2,j]       = rh_var[:,lat1:lat2,level,yearsInd[i,j]]
    u_tot_b[:,lat1:lat2,j]    = u_var[:,lat1:lat2,2,yearsInd[i,j]] - u_var[:,lat1:lat2,1,yearsInd[i,j]]
    v_tot_b[:,lat1:lat2,j]    = v_var[:,lat1:lat2,2,yearsInd[i,j]] - v_var[:,lat1:lat2,1,yearsInd[i,j]]
    VWS_tot_b[:,lat1:lat2,j]  = sqrt.(u_tot_b[:,lat1:lat2,j].^2 .+ v_tot_b[:,lat1:lat2,j].^2)
    #println("j value is: ",j,", and yearsInd is: ",yearsInd[i,j])
  end
  #println("ii value is: ",ii)
  # MPI values:
  mpi_test = mean(mpi_b, dims=3)
  mpi_yr[:,:,i] = mpi_test #should be the average of 6 fields...mpi_b[:,:,]
  # relative sst: rsst
  rsst_test = mean(sst_rel_b, dims=3)
  rsst_yr[:,:,i] = rsst_test #should be the average of 6 fields...mpi_b[:,:,]
  # RH values:
  rh_test = mean(rh_b, dims=3)
  rh_yr[:,:,i] = rh_test #should be the average of 6 fields...mpi_b[:,:,]
  # VWS values:
  vws_test = mean(VWS_tot_b, dims=3)
  vws_yr[:,:,i] = vws_test #should be the average of 6 fields...mpi_b[:,:,]
end

agrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
bgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
cgrid  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
cgrid2 = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
cgrid3 = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])
cgrid4 = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2])

for i in 1:dims[1]
  for j in 1:dims[2]
    agrid[i,j],bgrid[i,j] = find_best_fit(ri_events[:],mpi_yr[i,j,:])
    cgrid[i,j]            = cor(ri_events[:],mpi_yr[i,j,:])
    cgrid2[i,j]           = cor(ri_events[:],rsst_yr[i,j,:])
    cgrid3[i,j]           = cor(ri_events[:],rh_yr[i,j,:])
    cgrid4[i,j]           = cor(ri_events[:],vws_yr[i,j,:])
  end
end

levs = range(-2., 2., length = 21)
mpiReg = fig_anom_plot(agrid[:,:],lon,lat,"mpi linear regression with RI events",levs)
save("era5_mpi_reg.png", mpiReg, px_per_unit=6.0)

levs = range(-1., 1., length = 21)
vwsReg = fig_anom_plot(cgrid4[:,:],lon,lat,"vws linear regression with RI events",levs)
save("era5_vws_cor.png", vwsReg, px_per_unit=6.0)

levs = range(-1., 1., length = 21)
mpiCorr = fig_anom_plot(cgrid[:,:],lon,lat,"mpi correlation with RI events",levs)
save("era5_mpi_cor.png", mpiCorr, px_per_unit=6.0)

levs = range(-1., 1., length = 21)
sstCorr = fig_anom_plot(cgrid2[:,:],lon,lat,"rsst correlation with RI events",levs)
save("era5_sst_cor.png", sstCorr, px_per_unit=6.0)

levs = range(-1., 1., length = 21)
rhCorr = fig_anom_plot(cgrid3[:,:],lon,lat,"RH correlation with RI events",levs)
save("era5_rh_cor.png", rhCorr, px_per_unit=6.0)

levs = range(60., 90., length = 21)
# colormap = :vik, seems to work well for MPI
#mpiDiff = fig_anom_plot(mpi_b[:,:,1],lon,lat,"MPI (m/s), index 1",levs)
mpiDiff = fig_1_plot(mpi_b[:,:,1],lon,lat,"MPI (m/s), index 1",levs)
save("era5_mpi_1.png", mpiDiff, px_per_unit=6.0)

#mpiDiff = fig_1_plot(mpi_b[:,:,2],lon,lat,"MPI (m/s), index 1",levs)
#save("era5_mpi_2.png", mpiDiff, px_per_unit=6.0)
##
#mpiDiff = fig_1_plot(mpi_b[:,:,3],lon,lat,"MPI (m/s), index 1",levs)
#save("era5_mpi_3.png", mpiDiff, px_per_unit=6.0)
##
#mpiDiff = fig_1_plot(mpi_b[:,:,4],lon,lat,"MPI (m/s), index 1",levs)
#save("era5_mpi_4.png", mpiDiff, px_per_unit=6.0)
##
#mpiDiff = fig_1_plot(mpi_b[:,:,5],lon,lat,"MPI (m/s), index 1",levs)
#save("era5_mpi_5.png", mpiDiff, px_per_unit=6.0)
##
#mpiDiff = fig_1_plot(mpi_b[:,:,6],lon,lat,"MPI (m/s), index 1",levs)
#save("era5_mpi_6.png", mpiDiff, px_per_unit=6.0)
#
#mpiDiff = fig_1_plot(mpi_yr[:,:,1],lon,lat,"MPI (m/s), index 1",levs)
##mpiDiff = fig_1_plot(mpi_test[:,:,1],lon,lat,"MPI (m/s), index 1",levs)
#save("era5_mpi_yr_1.png", mpiDiff, px_per_unit=6.0)
#
#mpiDiff = fig_1_plot(mpi_yr[:,:,2],lon,lat,"MPI (m/s), index 1",levs)
##mpiDiff = fig_1_plot(mpi_test[:,:,1],lon,lat,"MPI (m/s), index 1",levs)
#save("era5_mpi_yr_2.png", mpiDiff, px_per_unit=6.0)

#mpiDiff = fig_1_plot(mpi_yr[:,:,12],lon,lat,"MPI (m/s), index 1",levs)
#save("era5_mpi_yr_12.png", mpiDiff, px_per_unit=6.0)










