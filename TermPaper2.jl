### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 52d41aa0-769c-11eb-2eb1-e78eeed9d5c8
begin
	using Plots
	using LsqFit
	using Statistics
	using PlutoUI
	md"_hidden code cell_"
end

# ╔═╡ c65eeff0-76b9-11eb-3543-3739a6484c24
html"""
<style>
.plutoui-toc.aside {
    top: 15%;
    bottom: 20%;
}</style>
<center><h1 style="font-size: 3.5rem">Term Paper 2<h1></center><br>
This work is a submission for the course Modelling Complex Systems(IDC-631) for the year 2020-21 Spring (Even) Semester by <b>Dhruva Sambrani, MS18163</b>
"""

# ╔═╡ d846746e-76bb-11eb-2938-adb0dc210f71
md"""# Abstract
This paper laregly recreates the paper[^soc] by Bak, Per and Tang, Chao and Wiesenfeld, Kurt on models that exhibit self-organized criticality. These models show certain structures repetitive structures and show distinctive power laws in certain properties. 

The specific model I have looked at is the Sandpile model, which has a relaxation time frequency having a log-log slope of between 1 to 2. This paper first introduces some working definitions and their equivalent definition in Julia[^julia] code.

First a pretty, but small simulation is shown which also verifies that the model is infact correctly translated. Post this, more rigourous analysis is done, and finally results and observations are mentioned
"""

# ╔═╡ a3a53a50-7830-11eb-380b-ed771520d060
md"## How to read the paper

This paper is made using the Pluto[^pluto] Package which provides an interface which can show both code and theory side by side, which will hopefully enhance the understanding of the work and also provide a reproducable implementation without having to go into the details of it. Those who wish to read the paper for the theory and results can easily skip over the code and still understand everything. Those who are attempting to recreate the work can use the code as a guide. However, in the interest of keeping the code aligned to the topic at hand, certain auxillary functions have been hidden. The content of these hidden cells can be found in the executable version of this notebook [here]()
"

# ╔═╡ f02e0ed0-784c-11eb-08b5-4746af933e1a
TableOfContents()

# ╔═╡ 2ed745de-76cf-11eb-101c-89f9d2ba170b
md"# Code Setup"

# ╔═╡ 1f57fcde-76cf-11eb-23c1-f5895772e792
md"## Packages used
The Plots.jl[^plots] package is used for plotting and the Statistics[^stats] Package for some statistical analysis. LsqFit.jl[^lsqfit] is used to fit the log-log plots to nice straight lines. PlutoUI.jl[^plutoui] is used for some UI elements such as the ToC
"

# ╔═╡ b239c390-76cf-11eb-1c95-4bb63c153f87
md"## Define Structure and Functions

A Sandpile grid is nothing but an Integer Array and doesn't need any special definitions. 
We just need to define the functions that act on a 2 dimensional grid which can simulate the rules of a Sandpile model
"

# ╔═╡ 1f872dbe-76d0-11eb-055c-cd0487ce223a
md"### The Sandpile model
In the most abstract sense, a Sandpile is just an array of integers which evolve in a specific way. However, more intuitively, it can be thought of as a distributions of equally sized grains of sand on a discrete grid. Then, the array of integers represent the height _or_ the slope of the pile in that location. Note that in the following report, we think of the values as the _slope_ of the pile. However, 
!!! warning \"The use of Height\"
    We use **height** to simply refer to the value of the cell, not the actual height of the pile. This poor choice of words is a common tactic to confuse readers so that they don't read the paper too carefully and find mistakes.

### Method of evolution

The model evolves with the following algorithm -

$$z_{i,j}(t+1) = z_{i,j}(t) - h_c(z_{i,j}(t)>h_c) + 
\text{count}(n\rightarrow h_c>n,\ \text{neighborhood})$$

Basically, we define a **neighborhood** of $h_c$ cells. If the height of the sand in a cell is greater than $h_c$, it is in a **critical** state, $h_c$ grains **topple** into the neighboring cells equally. Hence, we can simply subtract $h_c$ and count the number of **critical** neighbors to calculate the value of the cell in the next time step. Note that the sand is conserved in this form of evolution.

!!! info \"Relevant Quote\"
    Waste not the smallest thing for grains of sand make mountains

    **\- Eric Knight**
"

# ╔═╡ 5221fe40-7784-11eb-2d14-b7ca964ba546
md"""!!! tip "Our neighborhood and $h_c$" 
    In this work, we use the VonNeumann Neighborhood in 2D and radius 1. Hence, the neighbors are $(-1,0), (+1,0), (0,-1), (0,+1)$ and $h_c$ = 4
"""

# ╔═╡ b0797e20-769c-11eb-2711-7587647ce0c5
hc = 4

# ╔═╡ 7b2d6380-769c-11eb-12f6-35402dc5e96e
@inline iscritical(i::Int, hc::Int) = i >= hc

# ╔═╡ b1c3cc10-76d6-11eb-38e1-75a2c56e0044
md"#### Boundary Conditions

We use a fixed boundary condition of 0. This is like the edge of a table, because all sand keeps falling and is lost forever.

!!! info \"Relevant Quote\"
    Summer Vacation is slipping through our hands like grains of sand!

    **\- Calvin from Calvin and Hobbes**
"

# ╔═╡ 75a73440-769c-11eb-2001-a5345c6e3732
begin
	@inline function getindexbc(grid, i)
	    if checkbounds(Bool, grid, i)
	        @inbounds return grid[i]
	    else
	        return 0
	    end
	end
	getindexbc(grid) = i -> getindexbc(grid)
end

# ╔═╡ 8806b9d0-769c-11eb-28d9-bdc5b3870377
function evolveat(
        grid::Array{Int,2}, 
        hood::Array{CartesianIndex{2},1}, 
        index::CartesianIndex, 
        hc::Int,
        topple::Array{Int},
		affected_grid::BitArray{2},
        step::Int
    )
    res = grid[index]
    if res >= hc
        res -= 4
        @inbounds topple[step] += 1 # we'll talk about this later
		affected_grid[index] = true # and this
    end
    c = count(h -> iscritical(getindexbc(grid, h + index), hc), hood)
    res = res + c
    return res
end

# ╔═╡ 8dddef40-769c-11eb-2656-91b79e7af1cf
function evolve!(
        grid::Array{Int,2}, 
        simgrid::Array{Int,2},
		affected_grid::BitArray{2},
        hood::Array{CartesianIndex{2},1},
        hc::Int,
        topple::Array{Int,1},
        step::Int
    )
    for index in CartesianIndices(grid)
        @inbounds simgrid[index] = evolveat(grid, hood, index, hc, topple, affected_grid, step)
    end
	return simgrid, grid
end

# ╔═╡ 1c132690-76d9-11eb-0edf-45dda119d525
md"## Sanity Check and Pretty Animations

Let's run a small simulation on a 10x10 grid to see if we've coding in everything correctly. We start with an empty Pile and add a grain in the center of the grid. We wait for the pile to relax, and we add another grain of sand. We do this until we add  100 grains of sand.
"

# ╔═╡ 7db36160-76b8-11eb-22ec-bd23ba1800d6
let
	# This code is only to make the animation. The exact method of simulation is described later
	function anim_ready_next!(grid::Array{Int,2}, hc::Int, step::Int, max_steps::Int)
	    if count(iscritical.(grid, hc)) == 0
	        if step == max_steps
	            return false, step
	        else
	            grid[5, 5] += 1
	            step += 1
	        end
	    end
	    return true, step
	end
	function anim_simulate(grid::Array{Int}, hc::Int, max_steps::Int, res = false)
	    hood = CartesianIndex.([(-1, 0), (1, 0), (0, -1), (0, 1)])
	    step = 0
	    st = true
	    simgrid = similar(grid)
		results = Array{Int,2}[]
	    while st
	        st, step = anim_ready_next!(grid, hc, step, max_steps)
			grid, simgrid = evolve!(grid, simgrid, falses(size(grid)), hood, hc, [0], 0)
			push!(results, copy(grid))
	    end
	    return results
	end
	res = anim_simulate(zeros(Int, 10, 10), 4, 100);
	anim = @animate for i in res
		plot(i, st=:wireframe, fill=true, camera=(43,60), zlims=(0,4), clims=(0,4))
	end
	gif(anim, fps=5)
end

# ╔═╡ e0f4aab0-7797-11eb-0ecc-910d7b9ded81
md"_hidden code cell_"

# ╔═╡ f052c870-76d9-11eb-1b82-bb53ba216dc8
md"### Observations
1. The pile builds up to a height of 4 and then topples, as expected
2. The pile falls off at the edge, as expected
3. Some topplings are quick where as some topplings set off another topple and so on. However, these _chained_ topples are not as common as the single or double topples. We'll call an entire chain of topplings until a sand grain is added as a **avalanche**
4. This looks NOTHING like an actual sandpile. What were the authors thinking?\

Regarding point 4, this is because we are keeping a track of the _slope_, not the _height_ of the Sandpile. And yet we called it the _height_. I warned you, this is super confusing.
"

# ╔═╡ 0e96d220-76dc-11eb-2b86-f5708fe4803f
md"# Power law analysis

According to the paper mentioned before, there are certain properties which show a log-log frequency plot, which means that

$$\log(F_x) = m\cdot\log(x)+c$$

for certain properties of the model, where x is the value of the property and F(x) is the frequency of the value. These properties are -

1. **Number of sites affected** during an avalanche
2. **Total number of topples** induced during the avalanche
3. **Relaxation time** of the avalanche

Lets elucidate this further. Given a random condition, a single perturbation may lead to the pertubed cell to topple. 
This may, or may not, set off an avalanche which affects neighboring cells, has multiple topples, and takes a while to settle down. 
If we count these properties, and plot them on a log-log plot, it was shown that they form a straight line.
"

# ╔═╡ 08191482-7781-11eb-37bb-1fb9abbceee0
md"""!!! tip "Short note on Implementation"
    The properties are measured in the following fashion
    - 
        - Toggle count - The $i^{th}$ element of an array is updated everytime a topple is simulated in the $i^{th}$ perturbation.
        - Relaxation Time - The $i^{th}$ element of an array is updated after every `evolve` cycle for the $i^{th}$ perturbation.
        - Number of sites affected - The $(i,j)$ element in a grid of `false`s is set to `true` when the $(i,j)$ cell topples in a perturbation (not `evolve`) cycle. At the end of the $k^{th}$ perturbation, the total number of `true`s are counted and stored in the $k^{th}$ element of an array.
"""

# ╔═╡ 9484347e-769c-11eb-19c7-0961b3cdf2ef
function ready_next!(
		grid::Array{Int,2}, 
		affected_grid::BitArray{2},
		hc::Int,
		affect_count::Array{Int, 1},
		relax_time::Array{Int, 1},
		step::Int, max_steps::Int
	)
    if count(iscritical.(grid, hc)) == 0
		affect_count[step] = count(affected_grid)
		affected_grid[:] .= false
        if step == max_steps
            return false, step
        else
			# perturb
            xrand, yrand = rand.(Base.OneTo.(size(grid)))
            grid[xrand, yrand] += 1
            step += 1
        end
	else
		relax_time[step] += 1
	end
    return true, step
end

# ╔═╡ 9b1c21e0-769c-11eb-162e-35dbc2aba3f2
function simulate(grid::Array{Int}, hc::Int, max_steps::Int)
    hood = CartesianIndex.([(-1, 0), (1, 0), (0, -1), (0, 1)])
    topple_count = fill(0, max_steps)
	relax_time = fill(0, max_steps)
	affect_count = fill(0, max_steps)
    step = 1
    st = true
    simgrid = similar(grid)
	affected_grid = falses(size(grid))
    while st
        st, step = ready_next!(grid, affected_grid, hc, affect_count, relax_time, step, max_steps)
        grid, simgrid = evolve!(grid, simgrid, affected_grid, hood, hc, topple_count, step)
    end
    return topple_count, relax_time, affect_count
end

# ╔═╡ 36049990-778c-11eb-3007-4dc0bb6b2615
md"""## Method of analysis

Once we've simulated and collect the data, the analysis is done in the following way- 
1. Create a frequency table of `val => freq(val)` for each of the properties. 
2. Take a log for both values and frequencies.
3. Fit the data to $\text{freq} = m\cdot\text{val}+c$
4. Smooth the data and plot it against the fit

!!! warning "Things to keep in mind while fitting"
    The data does not follow a linear curve completely. The linearity is lost for reasons listed below. Hence, we should fit the data to only the first, linear part of the data, which usually is around 60%
"""

# ╔═╡ 668d2040-769d-11eb-3943-51fb559fded3
function freqtable(k::Array{T}) where T
    d = Dict{T,Int}()
    for i ∈ k
        if i != 0
            d[i] = get(d, i, 0) + 1
        end
    end
    return d
end

# ╔═╡ 6a004e02-769d-11eb-01c4-e5d635d9aa91
function getlogdata(topple_freq)
    x = log.(keys(topple_freq))
    perm = sortperm(x)
    y = log.(values(topple_freq))[perm]
    x = x[perm]
    return x, y
end

# ╔═╡ 6a009c20-769d-11eb-2d53-534435f34669
function smooth(xs, ys, length)
    binx = range(xs[1], stop=xs[end], length=length)
    _step = step(binx)
    biny = [
        mean(ys[findall(x -> val <= x < (val + _step), xs)])
        for val in binx
    ]
    binx = (binx .+ _step)[1:end - 1]
    return binx, biny
end

# ╔═╡ 4fd4eef0-769d-11eb-023e-21941d0cb482
function loglogfit(property_data)
	x, y = getlogdata(freqtable(property_data))
	smoothx, smoothy = smooth(x, y, 100)
	plot(smoothx,smoothy, st=:scatter)
	lnsp = floor(Int, 0.4*length(x))
	@. model(x, p) = p[1] * x + p[2]
	fit = curve_fit(model, x[1:lnsp], y[1:lnsp], [0.,0.])
	xp = 0:0.01:maximum(x)
	yp = model(xp, fit.param)
	return plot!(xp, yp, width=2), fit
end

# ╔═╡ 361d0af0-7789-11eb-29ca-d768d81ccc39
function stringifyfit(fit)
	val = fit.param[1]
	err = stderror(fit)[1]/val * 100
	"$(round(val, sigdigits=4))±$(round(err, sigdigits=4))%"
end

# ╔═╡ ef5f7640-7772-11eb-3d9e-4d7a092011df
function makeresults(all_data, gs, samples)
	tpl, rel, aff = loglogfit.(all_data)
	title!(tpl[1], 
		"Number of Topples\nm = $(stringifyfit(tpl[2]))",
		titlefontsize=10)
	title!(rel[1],
		"Relaxation Time\nm = $(stringifyfit(rel[2]))",
		titlefontsize=10)
	title!(aff[1],
		"Sites affected\nm = $(stringifyfit(aff[2]))",
		titlefontsize=10)
	plot(tpl[1], rel[1], aff[1], layout=(1,3), size=(3,1).*225)
	md" ## Grid size = $(gs), Samples = $(samples) $(current())"
end

# ╔═╡ 6224d450-7772-11eb-0efd-a96d742fa1b3
data_50_1_5 = simulate(rand(1:hc - 1, 50, 50), 4, 100000);

# ╔═╡ 0c34ab90-76a2-11eb-2612-071caa9a3824
makeresults(data_50_1_5, 50, "1e5")

# ╔═╡ 934a6d60-778b-11eb-2e1c-0d125b4b09df
data_50_1_4 = simulate(rand(1:hc - 1, 50, 50), 4, 10000);

# ╔═╡ bcaba3b0-76b7-11eb-1c8a-ff1e69e0e16e
makeresults(data_50_1_4, 50, "1e4")

# ╔═╡ abd6c270-778b-11eb-3d22-41ba454ab170
data_50_5_4 = simulate(rand(1:hc - 1, 50, 50), 4, 50000);

# ╔═╡ eca0f07e-7788-11eb-174d-fb38d141adb7
makeresults(data_50_5_4, 50, "5e4")

# ╔═╡ f4b73920-778b-11eb-3987-13806219df5e
data_20_5_4 = simulate(rand(1:hc - 1, 20, 20), 4, 50000);

# ╔═╡ 24b4ac10-7788-11eb-2b05-d3a7f4244703
makeresults(data_20_5_4, 20, "5e4")

# ╔═╡ 6f953350-779a-11eb-29f1-ab01cd638438
data_80_5_4 = simulate(rand(1:hc - 1, 80, 80), 4, 50000);

# ╔═╡ 5ab6ae90-779b-11eb-1727-9d2cda89e743
makeresults(data_80_5_4, 80, "5e4")

# ╔═╡ 795f8790-782c-11eb-2143-b9d7220990d3
begin
	f(dat) = loglogfit.(dat) .|> x->stringifyfit(x[2])
	mktable(i) = "| Sample Size | 20x20                | 50x50                | 80x80                |
| ----------- | -------------------- | -------------------- | -------------------- |
| 1e4         | No Data              | $(f(data_50_1_4)[i]) | No Data              |
| 5e4         | $(f(data_20_5_4)[i]) | $(f(data_50_5_4)[i]) | $(f(data_80_5_4)[i]) |
| 1e5         | No Data              | $(f(data_50_1_5)[i]) | No Data              |"
end

# ╔═╡ 21db9b60-784c-11eb-1a10-29a5f75821a4
md"## Result Summary"

# ╔═╡ 93a7faa0-7823-11eb-36e1-512899012714
Markdown.parse(
"### Topple Count
$(mktable(1))"
)

# ╔═╡ 7149d3ce-782c-11eb-2074-055c558ee198
Markdown.parse(
"### Relaxation time
$(mktable(2))"
)

# ╔═╡ d2393230-782c-11eb-1d53-dfaddfeb0ed2
Markdown.parse(
"### Affected grids
$(mktable(3))"
)

# ╔═╡ a6502560-7788-11eb-150a-cdf40bd5cbbd
md"## Observations

1. All the values of slope lie within the range $[-1.1, -0.95]$
2. As the sampling size increases, the error reduces
3. As the grid size increases, the error reduces
4. The linear domain 
    - increases with sample size until a point, after which it doesn't increase further
    - increases with grid size until a point, after which it doesn't increase further
"

# ╔═╡ 64b89820-782e-11eb-27df-c136babb8926
md"""# References
[^soc]: Self-organized criticality. Bak, Per and Tang, Chao and Wiesenfeld, Kurt. (Jul, 1988) Phys. Rev. A, 38: 364--374. doi: [10.1103/PhysRevA.38.364](https://link.aps.org/doi/10.1103/PhysRevA.38.364)
[^julia]: Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah. (2017) SIAM Review, 59: 65–98. doi: [10.1137/141000671](https://dx.doi.org/10.1137/141000671)
[^plots]: Plots - powerful convenience for visualization in Julia. Thomas Breloff. [docs.juliaplots.org](docs.juliaplots.org)
[^pluto]: Pluto - Reactive Notebooks for Julia. Fons van der Plas and Mikołaj Bochenski. [https://github.com/fonsp/Pluto.jl](https://github.com/fonsp/Pluto.jl)
[^plutoui]: PlutoUI. Fons van der Plas. [https://github.com/fonsp/PlutoUI.jl](https://github.com/fonsp/PlutoUI.jl)
[^stats]: Statistics - Julia stdlib for Statistics. [https://docs.julialang.org/en/v1/stdlib/Statistics/](https://docs.julialang.org/en/v1/stdlib/Statistics/)
"""

# ╔═╡ Cell order:
# ╟─c65eeff0-76b9-11eb-3543-3739a6484c24
# ╟─d846746e-76bb-11eb-2938-adb0dc210f71
# ╟─a3a53a50-7830-11eb-380b-ed771520d060
# ╟─f02e0ed0-784c-11eb-08b5-4746af933e1a
# ╟─2ed745de-76cf-11eb-101c-89f9d2ba170b
# ╟─1f57fcde-76cf-11eb-23c1-f5895772e792
# ╟─52d41aa0-769c-11eb-2eb1-e78eeed9d5c8
# ╟─b239c390-76cf-11eb-1c95-4bb63c153f87
# ╟─1f872dbe-76d0-11eb-055c-cd0487ce223a
# ╟─5221fe40-7784-11eb-2d14-b7ca964ba546
# ╠═b0797e20-769c-11eb-2711-7587647ce0c5
# ╠═7b2d6380-769c-11eb-12f6-35402dc5e96e
# ╠═8806b9d0-769c-11eb-28d9-bdc5b3870377
# ╟─b1c3cc10-76d6-11eb-38e1-75a2c56e0044
# ╠═75a73440-769c-11eb-2001-a5345c6e3732
# ╠═8dddef40-769c-11eb-2656-91b79e7af1cf
# ╟─1c132690-76d9-11eb-0edf-45dda119d525
# ╟─7db36160-76b8-11eb-22ec-bd23ba1800d6
# ╟─e0f4aab0-7797-11eb-0ecc-910d7b9ded81
# ╟─f052c870-76d9-11eb-1b82-bb53ba216dc8
# ╟─0e96d220-76dc-11eb-2b86-f5708fe4803f
# ╟─08191482-7781-11eb-37bb-1fb9abbceee0
# ╠═9484347e-769c-11eb-19c7-0961b3cdf2ef
# ╠═9b1c21e0-769c-11eb-162e-35dbc2aba3f2
# ╟─36049990-778c-11eb-3007-4dc0bb6b2615
# ╟─668d2040-769d-11eb-3943-51fb559fded3
# ╟─6a004e02-769d-11eb-01c4-e5d635d9aa91
# ╟─6a009c20-769d-11eb-2d53-534435f34669
# ╟─4fd4eef0-769d-11eb-023e-21941d0cb482
# ╟─361d0af0-7789-11eb-29ca-d768d81ccc39
# ╠═ef5f7640-7772-11eb-3d9e-4d7a092011df
# ╠═6224d450-7772-11eb-0efd-a96d742fa1b3
# ╟─0c34ab90-76a2-11eb-2612-071caa9a3824
# ╠═934a6d60-778b-11eb-2e1c-0d125b4b09df
# ╟─bcaba3b0-76b7-11eb-1c8a-ff1e69e0e16e
# ╠═abd6c270-778b-11eb-3d22-41ba454ab170
# ╟─eca0f07e-7788-11eb-174d-fb38d141adb7
# ╠═f4b73920-778b-11eb-3987-13806219df5e
# ╟─24b4ac10-7788-11eb-2b05-d3a7f4244703
# ╠═6f953350-779a-11eb-29f1-ab01cd638438
# ╟─5ab6ae90-779b-11eb-1727-9d2cda89e743
# ╟─795f8790-782c-11eb-2143-b9d7220990d3
# ╟─21db9b60-784c-11eb-1a10-29a5f75821a4
# ╟─93a7faa0-7823-11eb-36e1-512899012714
# ╟─7149d3ce-782c-11eb-2074-055c558ee198
# ╟─d2393230-782c-11eb-1d53-dfaddfeb0ed2
# ╟─a6502560-7788-11eb-150a-cdf40bd5cbbd
# ╟─64b89820-782e-11eb-27df-c136babb8926
