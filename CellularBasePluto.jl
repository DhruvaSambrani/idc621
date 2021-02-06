### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ fe6fc4e0-688c-11eb-3cba-fbadc8c45088
using Plots

# ╔═╡ b3150b20-68a7-11eb-0c1c-1389cd8227ec
html"""<center><h1 style="font-size: 3.5rem">CellularBase docs<h1></center><br>"""

# ╔═╡ ec8999c0-68a7-11eb-02dd-dfb6f5672595
md"This is the documentation of [CellularBase.jl](https://github.com/DhruvaSambrani/idc621/blob/main/CellularBase.jl)"

# ╔═╡ 37e3fad0-688c-11eb-0760-991b2ada6775
md"# Boundary Conditions

Define the kinds of Boundary Conditions

- Periodic
- FixedMin - Fixed BC with all elements of the min of the `possible_states`
- FixedMax - Fixed BC with all elements of the max of the `possible_states`
- Clamp - Nearest cell within the region of interest

"

# ╔═╡ 29996140-688c-11eb-1317-6d36dc7fe265
@enum BoundaryCondition Periodic FixedMin FixedMax Clamp

# ╔═╡ a222c570-688c-11eb-0a24-bf1da8cab724
md"# Neighborhoods
Every neighborhood has a `radius` and a set of `CartesianIndex`es associated with it. The most general way of evaluating these is by referencing the appropriate property.
"

# ╔═╡ 2ddfe490-688c-11eb-11fe-7fe9cf9cc425
begin
	# Neighborhoods
	abstract type AbstractNeighborhood end
	radius(neighborhood::AbstractNeighborhood) = neighborhood.radius
	cartesians(neighborhood::AbstractNeighborhood) = neighborhood.cartesians
end

# ╔═╡ eabf0730-688c-11eb-25cd-0fdf3696c423
md"## VonNeumann Neighborhood

The VonNeumann Neighborhood is the simplest neighborhood which contains $k$ neighbors adjacent to it.

It is defined for multiple dimensions.

A 2D representation of radius(=$k$) 3 is shown below
"

# ╔═╡ a1003e70-688c-11eb-12c6-193712a896ff
struct VonNeumannNeighborhood <: AbstractNeighborhood
	cartesians::Array{CartesianIndex}
	radius::Int
	dimensions::Int
	function VonNeumannNeighborhood(radius::Int, dimensions::Int)::VonNeumannNeighborhood
		arr = CartesianIndex[CartesianIndex(Tuple(zeros(Int, dimensions)))]
		for i in 1:radius
			for j in 1:dimensions
				zs = zeros(Int, dimensions)
				zs[j] = i
				push!(arr, CartesianIndex(Tuple(zs)))
				zs[j] = -i
				push!(arr, CartesianIndex(Tuple(zs)))
			end
		end
		new(arr, radius, dimensions)
	end
end

# ╔═╡ af3a67d0-688d-11eb-1d4d-bff9ccfbd737
begin
	vline([-4, 4], color = :black)
	hline!([-4, 4], color = :black)
	scatter!(Tuple.(cartesians(VonNeumannNeighborhood(3, 2))), grid=false, legend=false)
	scatter!((0,0), size=(240, 160), axis = false, ticks = false, lims = (-4,4))
end

# ╔═╡ 4241b7e0-688e-11eb-3c30-f7c0f1730888
md"## Moore Neighborhood

The Moore Neighborhood is a simple neighborhood which contains all neighbours in the $ n$ d-cube centered at the cell and of side length $2k+1$

It is defined for multiple dimensions.

A 2D representation of radius(= $k$) 3 is shown below
"

# ╔═╡ 3752f380-688e-11eb-31bd-699522b02cbd
struct MooreNeighborhood <: AbstractNeighborhood
	cartesians::Array{CartesianIndex}
	radius::Int
	dimensions::Int
	function MooreNeighborhood(radius::Int, dimensions::Int)::MooreNeighborhood
		new(
			CartesianIndices(Tuple(-radius:radius for _ in 1:dimensions)) |>
			x -> reshape(x, :),
			radius, dimensions
		)
	end
end

# ╔═╡ d069c712-688e-11eb-38c1-47e1dec2e2aa
begin
	vline([-4, 4], color = :black)
	hline!([-4, 4], color = :black)
	scatter!(Tuple.(cartesians(MooreNeighborhood(3, 2))), grid=false, legend=false)
	scatter!((0,0), size=(240, 160), axis = false, ticks = false, lims = (-4,4))
end

# ╔═╡ 61003e30-688f-11eb-1732-3dc61a8abad2
md"# Grids
These are the actual cellular automata which are the aim of this package.

`AbstractGrid` is concretized to a working grid which has the real working mechanism
"

# ╔═╡ 43d5cff0-688f-11eb-0cbe-d32c52e694cf
abstract type AbstractGrid{T} end

# ╔═╡ a86140ce-688f-11eb-21eb-c3cd76fcdcd7
md"## `getindex`
`Base.getindex` is specialized to AbstractGrids to handle Boundary conditions
"

# ╔═╡ e27b1022-688f-11eb-0aee-85fb8c9d3834
md"""## State functions

These are functions related to the _state_ of the grid.

!!! warning "Optimization Possibility"
    `newstate` and `state!` can often be specialized and optimized to be more efficient and hence throw warnings
"""

# ╔═╡ f6fb04b0-688f-11eb-16a7-659171a260de
begin
	function state(grid::AbstractGrid{T})::Array{T} where {T} 
		grid.state
	end
	
	function state!(grid::AbstractGrid{T}, newstate)::Nothing  where {T}
		@warn "Possible unspecialized call"
		grid.state = newstate
		return nothing
	end
	
	function newstate(grid::AbstractGrid{T})::Array{T} where {T} 
		@warn "Possible unspecialized call"
		similar(state(grid))
	end
end;

# ╔═╡ 1f21cbe0-6890-11eb-309b-a7512538a923
md"## Accessor functions"

# ╔═╡ feac0330-688f-11eb-03a7-f7039760a343
begin
	boundary_conditions(grid::AbstractGrid)::BoundaryCondition = grid.bc
	neighborhood(grid::AbstractGrid)::Array{CartesianIndex} = cartesians(grid.neighborhood)
	function possible_states(grid::AbstractGrid{T})::Array{T} where {T}
		grid.possible_states
	end
	function neighbors(grid::AbstractGrid{T}, i::CartesianIndex)::Array{T} where {T}
		neighborindex = neighborhood(grid) .|> x -> x + i
		grid[neighborindex]
	end
end;

# ╔═╡ 5b49fe40-688f-11eb-333c-c722511ca6e0
function Base.getindex(grid::AbstractGrid{T}, index::CartesianIndex)::T where {T}
    try
        state(grid)[index]
    catch
        if boundary_conditions(grid) == Periodic
            grid[CartesianIndex(mod1.(Tuple(index), size(grid)))] # mod-1 arithmetic
        elseif boundary_conditions(grid) == FixedMax
            maximum(possible_states(grid))
        elseif boundary_conditions(grid) == FixedMin
            minimum(possible_states(grid))
        else
            grid[
                CartesianIndex(
                    Tuple(
                        clamp(Tuple(index)[i], 1, size(grid)[i]) 
                        for i in 1:ndims(grid)
                    )
                )
            ]
        end
    end
end

# ╔═╡ c6e11070-6890-11eb-3f5a-051919a603b2
md"## Useful Helpers"

# ╔═╡ ea9a4450-6890-11eb-334e-0d4d584222db
begin
	Base.ndims(grid::AbstractGrid) = ndims(state(grid))
	Base.size(grid::AbstractGrid) = size(state(grid))
	function Base.getindex(grid::AbstractGrid{T}, indexs::Union{Array{Int},Array{CartesianIndex{N}}})::Array{T} where {T,N}
		indexs .|> i -> grid[i]
	end
	
	function Base.getindex(grid::AbstractGrid{T}, is::Int...)::T where T
		grid[CartesianIndex(is)]
	end
	Base.CartesianIndices(grid::AbstractGrid)::CartesianIndices = CartesianIndices(state(grid))
end

# ╔═╡ 0c8a9f10-6891-11eb-044c-abe0b2b66629
md"# Evolutions and Simulations
This section defines the most abstract functions which evolve and simulate a Cellular Automaton
"

# ╔═╡ 4a177510-6891-11eb-12ca-35c2f8dad57d
md"""## Evolutions

`evolve!` is a mutator function that evolves the grid by one step. A new state is created by calling `newstate`. Each index is iterated through, and the values of the neighbours are sent as parameters to the `grid`. The grid then evaluates the next state of the present index, and hence the evolution rules must be coded into the `grid` as a functional object

!!! danger "Grid is a Function too!"
    Every AbstractGrid is also assumed to be a function-like object and the code will throw an error if it is not defined to be as such

The new state is then set back into the grid.

`tabular_evolution!` is a helper function that can be used for rules which can be encoded into a table, like Wolfram's CA.

!!! info "tabular_evolution!"
    This function can be used in a specialized `evolve!`. Send the table via the `kwargs` to `evolve!` and then unpack to `tabular_evolution!`. See [caWolfram.jl](https://github.com/DhruvaSambrani/idc621/blob/main/caWolfram.jl) for an implementation.

"""

# ╔═╡ 2ea1bbb0-6891-11eb-3c44-2706f3faa991
begin
	function evolve!(grid::AbstractGrid{T}; kwargs...)::Nothing where {T}
	    _newstate = newstate(grid)
	    Threads.@threads for i in CartesianIndices(grid)
	        _newstate[i] = grid(neighbors(grid, i); kwargs...)
		end
		state!(grid, _newstate)
		return nothing
	end
	
	function tabular_evolution!(grid::AbstractGrid, table::Dict)::Nothing
	    newstate = newstate(grid)
	    for i in CartesianIndices(grid)
	        newstate[i] = table[grid(neighbors(grid, i))]
		end
		state!(grid, newstate)
		return nothing
	end
end

# ╔═╡ bc977ad0-6892-11eb-065f-8561963d653e
md"""## Simulation
The outer most function that will almost always be called _unspecialized_ by the user.
The `simulate!` function takes a grid and a number of steps to simulate for. The `store_results` parameter allows users to not store results after every `evolve!` call. `postrunhook` is called after every step with parameters as `(step, grid, kwargs...)`. All `kwargs` are forwarded to `evolve!` and `postrunhook`.

!!! tips "postrunhook Special Return Values"
    While `postrunhook` can mutate the grid on specific steps, it can also act as a messenger to `simulate!` to change its behavior. See below for available messages

List of special messages-

1. `:shortcurcuit` ⟹ do not call `evolve!` on the grid in the next step. This can be useful if you want to stop evolution for some steps, or if you know that the grid will not change anymore, such as when it has reached a stable state.
2. `:interrupt` ⟹ return immediately.

!!! info "Implementation Example"
    See [Term Paper 1#Equilibrium infection distribution](https://dhruvasambrani.github.io/idc621/termpaper1) for example

"""

# ╔═╡ e450da2e-6892-11eb-21e4-c3e5d6141a6d
function simulate!(grid::AbstractGrid{T}, steps::Int; store_results=true, postrunhook=nothing, kwargs...)::Array{Array{T}}  where {T}
	if store_results results = Array{T}[] end
	push!(results, copy(state(grid)))
	if !isnothing(postrunhook) 
		postrunhook(grid, step; kwargs...)
	end
	ss = false
	for step in 1:steps
		if (!ss) evolve!(grid; kwargs...) end
		if store_results push!(results, copy(state(grid))) end
		if !isnothing(postrunhook)
			message = postrunhook(grid, step; kwargs...)
			ss = message == :shortcurcuit
			if message == :interrupt
				break
			end
		end
	end
	if store_results results else nothing end
end

# ╔═╡ Cell order:
# ╟─b3150b20-68a7-11eb-0c1c-1389cd8227ec
# ╟─ec8999c0-68a7-11eb-02dd-dfb6f5672595
# ╠═fe6fc4e0-688c-11eb-3cba-fbadc8c45088
# ╟─37e3fad0-688c-11eb-0760-991b2ada6775
# ╠═29996140-688c-11eb-1317-6d36dc7fe265
# ╟─a222c570-688c-11eb-0a24-bf1da8cab724
# ╠═2ddfe490-688c-11eb-11fe-7fe9cf9cc425
# ╟─eabf0730-688c-11eb-25cd-0fdf3696c423
# ╟─af3a67d0-688d-11eb-1d4d-bff9ccfbd737
# ╠═a1003e70-688c-11eb-12c6-193712a896ff
# ╟─4241b7e0-688e-11eb-3c30-f7c0f1730888
# ╟─d069c712-688e-11eb-38c1-47e1dec2e2aa
# ╠═3752f380-688e-11eb-31bd-699522b02cbd
# ╟─61003e30-688f-11eb-1732-3dc61a8abad2
# ╠═43d5cff0-688f-11eb-0cbe-d32c52e694cf
# ╟─a86140ce-688f-11eb-21eb-c3cd76fcdcd7
# ╠═5b49fe40-688f-11eb-333c-c722511ca6e0
# ╟─e27b1022-688f-11eb-0aee-85fb8c9d3834
# ╠═f6fb04b0-688f-11eb-16a7-659171a260de
# ╟─1f21cbe0-6890-11eb-309b-a7512538a923
# ╠═feac0330-688f-11eb-03a7-f7039760a343
# ╟─c6e11070-6890-11eb-3f5a-051919a603b2
# ╠═ea9a4450-6890-11eb-334e-0d4d584222db
# ╟─0c8a9f10-6891-11eb-044c-abe0b2b66629
# ╟─4a177510-6891-11eb-12ca-35c2f8dad57d
# ╠═2ea1bbb0-6891-11eb-3c44-2706f3faa991
# ╟─bc977ad0-6892-11eb-065f-8561963d653e
# ╠═e450da2e-6892-11eb-21e4-c3e5d6141a6d
