#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test functions, mostly taken from juliadatascience.pdf
#
# run with: 
# /Users/silvers/.juliaup/bin/julia functest.jl
#
# to run from within a juila REPL: 
# include("/Users/silvers/code/juliaCode/functest.jl")
#
# once functest.jl has been run in the REPL, then these functions can be
# accessed within the REPL, e.g.: 
# add_numbers(2, 3) --> 5
# figure_canvas()   --> creates and opens a png figure.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using DataFrames
using CSV
using CairoMakie
using Statistics
using NCDatasets


function add_numbers(x, y)
    return x + y
end

#add_numbers(17, 29)

function grades_array()
    name = ["Bob", "Sally", "Alice", "Hank"]
    age  = [17, 18, 20, 19]
    grade_2020 = [5.0, 1.0, 8.5, 4.0]
    (; name, age, grade_2020)
end

function grades_2020()
    name2 = ["Sally", "Bob", "Alice", "Hank"]
    grade_20202 = [1, 5, 8.5, 4]
    DataFrame(; name2, grade_20202)
end
grades_2020()


function figure_canvas()
    fig = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(600,400),
        )
        ax = Axis(fig[1,1];
            xlabel="x",
            ylabel="y",
            title="Title",
        ) 
    ylims!(-1, 2)
    xlims!(0, 13)
    lines!(ax, 0.5:0.2:4pi, x -> cos(x)/x;
        color=:black,
        linewidth=2,
        linestyle=:dash,
        label = "cos(x)/x",
        )
    scatterlines!(ax, 0.5:0.2:4pi, x -> -cos(x)/x;
        color=:black,
        linewidth=2,
        linestyle=:dash,
        label = "-cos(x)/x",
        )
    axislegend("legend"; position=:rt)
    fig # does this send the plot to the png file?
end
#figure_canvas()

## code trying to work with DataFrames and ENSO indices...
#dfe = CSV.read(fileENSO, header = 4, delim="  ", footerskip = 4, DataFrame)
#nms = ["year", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
#dff = DataFrame(dfe, nms)
#equals_yr(year::String7) = year == " 2005"
#filter(:year => equals_yr, dff)

function loop()
    for i in 1:10
      println(i)
      equals_yr(year::String7) = year == " 2005"
    end
end

function testMult(x, y)
    return x * y
end

function multiple_scatters_and_lines()
    x = collect(0:10)
    cycle = Cycle([:color, :linestyle, :marker], covary=true) 
    set_theme!(Lines=(cycle=cycle,), Scatter=(cycle=cycle,)) 
    fig = Figure(size=(600, 400), font="CMU Serif")
    ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"f(x,a)")
    for i in x
        lines!(ax, x, i .* x; label=L"%$i x")
        scatter!(ax, x, i .* x; markersize=13, strokewidth=0.25, 
            label=L"%$i x")
    end
    axislegend(L"f(x)"; merge=true, position=:lt, nbanks=2, labelsize=14)
    text!(L"f(x,a) = ax", position=(4, 80))
    set_theme!() # reset to default theme
    fig
end

function set_colors_and_cycle()
#  see section 6.7.1 of juliadatascience.pdf for this function.  
#    # Epicycloid lines
#    # how do I write greek letters in vim? 
#    x(r, k, \theta) = 
#    x(r, k, θ) = r∗(k .+ 1.0) .∗cos.(θ) .− r∗cos.((k .+ 1.0) .∗θ) 
#    y(r, k, θ) = r ∗ (k .+ 1.0) .∗ sin.(θ) .− r ∗ sin.((k .+ 1.0) .∗ θ) 
#    θ = range(0, 6.2π, 1000)
#    
#    axis = (; xlabel=L"x(\theta)", ylabel=L"y(\theta)",
#        title="Epicycloid", aspect=DataAspect())
#    figure = (; size=(600, 400), font="CMU Serif")
#    
#    fig, ax, _ = lines(x(1, 1, θ), y(1, 1, θ); color="firebrick1", # string
#       label=L"1.0", axis=axis, figure=figure)
#    lines!(ax, x(4, 2, θ), y(4, 2, θ); color=:royalblue1, #symbol
#        label=L"2.0")
#    for k = 2.5:0.5:5.5
#        lines!(ax, x(2k, k, θ), y(2k, k, θ); label=latexstring("$(k)")) #cycle
#    end
#    Legend(fig[1, 2], ax, latexstring("k, r = 2k"), merge=true)
#    colsize!(fig.layout, 1, Aspect(1, 1.0))
#    fig
end

# to use: filter([:YEAR, :WIND_KT] => first_tc_filter, dfRI)
function first_tc_filter(year, wind)::Bool
    test1 = 1999 < year
    test2 = 80 < wind
    test1 && test2
end

function first_nina_filter(year, basin)::Bool
    # NH La Nina years: 1998, 1999, 2007, 2010, 2016, 2020
    #                   2021, 2022
    # NH El Nino years: 1991, 1994, 1997, 2002, 2004, 2009,
    #                   2015, 2023
    #test1 = 2005 == year
    #nope test1 = (2005 || 2019) == year
    test1 = (year == 1998) || (year == 1999) || (year == 2007) || (year == 2010) || (year == 2016) || (year == 2020) || (year == 2021) || (year == 2022)
    test2 = basin == "NA"
    test1 && test2
end

function first_nino_filter(year, basin)::Bool
    test1 = (year == 1991) || (year == 1994) || (year == 1997) || (year == 2002) || (year == 2004) || (year == 2009) || (year == 2015) || (year == 2023)
    test2 = basin == "NA"
    test1 && test2
end

function sec_nino_filter(year, basin, basString::String)::Bool
    test1 = (year == 1991) || (year == 1994) || (year == 1997) || (year == 2002) || (year == 2004) || (year == 2009) || (year == 2015) || (year == 2023)
    test2 = basin == basString
    test1 && test2
end

function thd_nino_filter(year)::Bool
    test1 = (year == 1991) || (year == 1994) || (year == 1997) || (year == 2002) || (year == 2004) || (year == 2009) || (year == 2015) || (year == 2023)
    #test2 = basin == "NA"
    #test1 && test2
end
function thd_nina_filter(year)::Bool
    test1 = (year == 1998) || (year == 1999) || (year == 2007) || (year == 2010) || (year == 2016) || (year == 2020) || (year == 2021) || (year == 2022)
end 

function make_hist(dfin,tit)
    fig = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:snow2,
        size=(400,400),
        )
    ax = Axis(fig[1,1];
        xlabel="Max Wind(kt) at Start of RI",
        ylabel="# of storms",
        title=tit,
        )
    hist!(ax, dfin, bins=35:5:100)
    fig
end

function cshift(lons, field, lon_0)
   shift = @. lons - lon_0 > 180
   nn = sum(shift)
   (circshift(lons - 360shift, nn), circshift(field, (nn, 0)))
end








