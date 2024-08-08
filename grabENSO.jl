
equals_yr(year::String7) = year == " 2005"

whichy = [1965, 1970, 1975, 1980, 1985]

for i in 1:4
   print(whichy[i]), print(" ")
   equals_yr(year::String7) = year == " 2005"
end

# grap a row...
boo = dfe[4, 2:13]

# I want boo to be a vector, composed of the jan:dec columns... 
# and then I want to append many vectors onto each other so that 
# I can end up with a long vector that includes many years one
# after another.  

# but boo is still not a vector, and I don't know how to append
# vectors.  

# if we could put this function in a loop that cycles through the
# rows of choice, and add up the result of each iteration, I think
# I might be happy.  
function grab_a_row(i::String7)
    j = parse(Int64, i)
    dfe[j, :]
end

c = [a; b]
