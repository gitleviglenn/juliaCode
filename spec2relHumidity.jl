#----------------------------------------------------------------------------
# calculate the relative humidity using the specific humidity and temperature.   
# this does not appear to be working.  i don't trust the fields that 
# are produced....
#
# levi silvers                                             June 2025
#----------------------------------------------------------------------------

function spec2rel(filein)
  # open file and read specific humidity data and temperature data
  # compute the relative humidity from the specific
  return relH
end

function testFun(num)
  result = 5 * num
  return result
end

function spec2relComp(path,filein,ind_nino,ind_nina)

  file2    = path*"MPI_CESM2-FV2_full_output.nc"
  
  scaleF   = 1000000
  
  data2 = NCDataset(file2)
  #vmax  = data2["vmax"]
  shum  = data2["q"] # specific humidity (144, 73, 19, 1032) (kg/kg)
  temp  = data2["t"] # air temperature (K)
  lev   = data2["plev"]
  
  #sst_var     = data2["shum"]
  dims        = size(sst_var)
  dims2       = size(ind_nino) 
  println("levels are: ",lev[:])
  println("level 4 is: ",lev[4])

  timesteps = dims2[1] 
  # change numfields to timesteps 
  temp_trH  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], timesteps);
  temp_trL  = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], timesteps);
  q_trH     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], timesteps);
  q_trL     = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], timesteps);
  num       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], timesteps);
  den       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], timesteps);
  svp       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], timesteps);
  rhH       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], timesteps);
  rhL       = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], timesteps);
  global relHcomp  = Array{Union{Missing, Float64}, 2}(undef, dims[1], dims[2]);
  
  println("size of temp is: ",size(temp))
  println("size of temp_tr is: ",size(temp_trH))
  for i in 1:timesteps
    temp_trH[:,:,i] = temp[:,:,4,ind_nino[i]] .- 273.15 # convert to C   
    q_trH[:,:,i]    = shum[:,:,4,ind_nino[i]]   
    temp_trL[:,:,i] = temp[:,:,4,ind_nina[i]] .- 273.15 # convert to C 
    q_trL[:,:,i]    = shum[:,:,4,ind_nina[i]]   
  end
  
  #println("temperature is (should be celcius): ",temp_trH[:,10,1])
  println("specific humidity is (should be kg/kg): ",q_trH[:,10,1])
  
  tempT = temp_trH
  qT    = q_trH
  num = exp.(34.494 .- 4924.99./(tempT .+ 237.1))
  den = (tempT .+ 105).^1.57
  svp = num./den
  # reuse the num and den arrays:
  num = qT 
  den = (1. .- qT).*svp
  # relative humidity from high ENSO values
  rhH  = num./den
  
  tempT = temp_trL
  qT    = q_trL
  num = exp.(34.494 .- (4924.99./(tempT .+ 237.1)))
  den = (tempT .+ 105).^1.57
  svp = num./den
  # reuse the num and den arrays:
  num = qT 
  den = (1. .- qT).*svp
  # relative humidity from low ENSO values
  rhL  = num./den
  
  # compute the composite relative humidity field
  rh_comp = rhH .- rhL
 
  # it seems like no error message is generated if the dims variable passed to mean() 
  # is nonsense....   why?    
  # these magnitudes don't make sense to me at all...
  #RelRH_comp_mn = scaleF .* mean(rh_comp, dims=4)
  println("dimensions of rhH and rh_comp: ",size(rhH)," ",size(rh_comp))
  println("number of timesteps: ",size(ind_nino))
  relHcomp = 1 * scaleF .* mean(rh_comp, dims=3)

  return relHcomp
end
