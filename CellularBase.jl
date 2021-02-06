module CellularBase
import Base: getindex, ndims, size, CartesianIndices, show
import Markdown: parse

export BoundaryCondition, Periodic, FixedMax, FixedMin, Clamp, boundary_conditions,
		AbstractNeighborhood, VonNeumannNeighborhood, MooreNeighborhood,
		radius, cartesians, neighborhood,
		AbstractGrid, state, state!, possible_states, newstate,
		tabular_evolution!, simulate!

@enum BoundaryCondition Periodic FixedMin FixedMax Clamp

# Neighborhoods
abstract type AbstractNeighborhood end
radius(neighborhood::AbstractNeighborhood) = neighborhood.radius
cartesians(neighborhood::AbstractNeighborhood) = neighborhood.cartesians

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

function Base.show(io::IO, ::MIME"text/html", neighborhood::AbstractNeighborhood)
	show(
		io,
		MIME("text/html") ,
		parse("**Type:** " * repr(typeof(neighborhood), context=io) *
		"\n\n **Radius:** " * repr(radius(neighborhood), context=io) *
		"\n\n **List of Neighbors:** " * repr(Tuple.(cartesians(neighborhood)), context=io))
	)
end

# Grids
abstract type AbstractGrid{T} end

function Base.getindex(grid::AbstractGrid{T}, index::CartesianIndex)::T where {T}
    try
        state(grid)[index]
    catch
        if boundary_conditions(grid) == Periodic
            grid[CartesianIndex(mod1.(Tuple(index), size(grid)))]
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

function Base.show(io::IO, ::MIME"text/html", grid::AbstractGrid)
	println(
		io,
		"Type: ", typeof(grid),
		"Possible States: ", possible_states(grid),
		"<br>Boundary Condition: ", boundary_conditions(grid),
		"<br>Neighborhood: ", neighborhood(grid),
		"<br>Size: ", size(grid)
	)
end

function state(grid::AbstractGrid{T})::Array{T} where {T} 
	grid.state
end

function state!(grid::AbstractGrid{T}, newstate)::Nothing  where {T}
	@warn "Possible unspecialized Call"
	grid.state = newstate
	return nothing
end
boundary_conditions(grid::AbstractGrid)::BoundaryCondition = grid.bc
neighborhood(grid::AbstractGrid)::Array{CartesianIndex} = cartesians(grid.neighborhood)
function possible_states(grid::AbstractGrid{T})::Array{T} where {T}
	grid.possible_states
end
function newstate(grid::AbstractGrid{T})::Array{T} where {T} 
	@warn "Possible unspecialized Call"
	similar(state(grid))
end
function neighbors(grid::AbstractGrid{T}, i::CartesianIndex)::Array{T} where {T}
    neighborindex = neighborhood(grid) .|> x -> x + i
	grid[neighborindex]
end

Base.ndims(grid::AbstractGrid) = ndims(state(grid))
Base.size(grid::AbstractGrid) = size(state(grid))
function Base.getindex(grid::AbstractGrid{T}, indexs::Union{Array{Int},Array{CartesianIndex{N}}})::Array{T} where {T,N}
	indexs .|> i -> grid[i]
end

function Base.getindex(grid::AbstractGrid{T}, is::Int...)::T where T
	grid[CartesianIndex(is)]
end
Base.CartesianIndices(grid::AbstractGrid)::CartesianIndices = CartesianIndices(state(grid))

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

end
