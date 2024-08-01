#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotENSO.jl
#
# - plot an ENSO index verses time
# - initially the nino3.4 index is used but this script should be 
# - generalized to plot ENSO according to mulitple indices
#
# to do: 
#     - incorporate additiona ENSO indices
#     - dramatically improve the code that concatinates the data, i 
#       clearly know nothing about julia at this time 
#
# levi silvers                                             august 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
using DataFrames
using CSV
using CairoMakie

# run with: 
# /Users/silvers/.juliaup/bin/julia script.jl

# data over 73 years
fileENSO = "/Users/silvers/data/enso/nina34.anom.data"    

dfe = CSV.read(fileENSO, header = 4, delim="  ", footerskip = 4, DataFrame)  

nms = ["year", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

dff = DataFrame(dfe, nms)

# filter may also work well for selecting rows...
dff.t = parse.(Float64, dff.jan)
janTest  = transpose(dff.t)     # produces a 1x73 element vector
febTest  = transpose(parse.(Float64, dff.feb)) # produces a 1x73 elem vect
marTest  = transpose(parse.(Float64, dff.mar)) # produces a 1x73 elem vect
aprTest  = transpose(parse.(Float64, dff.apr)) # produces a 1x73 elem vect
mayTest  = transpose(parse.(Float64, dff.may)) # produces a 1x73 elem vect
junTest  = transpose(parse.(Float64, dff.jun)) # produces a 1x73 elem vect
julTest  = transpose(parse.(Float64, dff.jul)) # produces a 1x73 elem vect
augTest  = transpose(parse.(Float64, dff.aug)) # produces a 1x73 elem vect
sepTest  = transpose(parse.(Float64, dff.sep)) # produces a 1x73 elem vect
octTest  = transpose(parse.(Float64, dff.oct)) # produces a 1x73 elem vect
novTest  = transpose(parse.(Float64, dff.nov)) # produces a 1x73 elem vect
decTest  = transpose(parse.(Float64, dff.dec)) # produces a 1x73 elem vect

#id = [1, 2, 3]
#tdf = DataFrame(; id, janTest, febTest, marTest)

fig = Figure()

ax = Axis(fig[1,1];
    xlabel="monthly mean values", 
    ylabel="anomaly",
    title="ENSO index over 12 years"
    )

index = 1
for i in 1:10
    println(i)
    println(janTest[index])
    println(janTest[i])
end

ensoTS = [janTest[index], febTest[index], marTest[index], aprTest[index], mayTest[index], junTest[index], julTest[index], augTest[index], sepTest[index], octTest[index], novTest[index], decTest[index]]


ensoTSa = Vector{Float64}(undef, 12)
someone = Vector{Float64}(undef, 100)

#i = 1
#for i in 1:10
#  global ensoTSa = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
#end

# there has to be a better way to do this with a simple loop
index=61
i = index+1
ensoTSa = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
i = index+2
ensoTSb = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
i = index+3
ensoTSc = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
i = index+4
ensoTSd = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
i = index+5
ensoTSe = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
i = index+6
ensoTSf = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
i = index+7
ensoTSg = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
i = index+8
ensoTSh = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
i = index+9
ensoTSi = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
i = index+10 
ensoTSj = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
i = index+11 
ensoTSk = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]
i = index+12 
ensoTSl = [janTest[i], febTest[i], marTest[i], aprTest[i], mayTest[i], junTest[i], julTest[i], augTest[i], sepTest[i], octTest[i], novTest[i], decTest[i]]

ensocat = cat(ensoTSa, ensoTSb, ensoTSc, ensoTSd, ensoTSe, ensoTSf, ensoTSg, ensoTSh, ensoTSi, ensoTSj, ensoTSk, ensoTSl, dims = 1)

println(ensoTSa)
describe(ensoTSa)
describe(someone)
describe(ensocat)

#lines!(ax, dff.t, linewidth = 3.0)
#lines!(ax, marTest[:], linewidth = 3.0)
#lines!(ax, octTest[:], linewidth = 3.0)
lines!(ax, ensocat[:], linewidth = 3.0)

#lines!(ax, dff[1, 2:13], linewidth = 3.0)

save("plotENSO.png",fig)

# careful, I think this will add to the dataframe dff!!

#dff.mjuly = parse.(Float64, dff.jul) # produces a 73 element vector
#julyTest  = transpose(dff.mjuly)     # produces a 1x73 element vector
# creates a 73 element vector...

# the above script will plot enso values for every value of january over the 
# 73 year period.  Add the defualt vectors that I get together would result 
# in plotting all the january's first, then all the february's, then all the 
# marches, etc....   

# what I really want is to plot the enso values for a given year, and then add
# the desired years up so that I can plot a time series of enso between 1980 and 
# 2024.  i think the transpose operation will all this.  






