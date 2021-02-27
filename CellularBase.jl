module CellularBase
import Base: getindex, ndims, size, CartesianIndices, show, eltype
import Markdown:parse

export BoundaryCondition, Periodic, FixedMax, FixedMin, Clamp, boundary_conditions,
		AbstractNeighborhood, VonNeumannNeighborhood, MooreNeighborhood,
		radius, cartesians, neighborhood,
		AbstractGrid, state, state!, possible_states, newstate,
		evolve!, tabular_evolution!, simulate!

abstract type BoundaryCondition end
struct Periodic <: BoundaryCondition end
struct FixedMin <: BoundaryCondition end
struct FixedMax <: BoundaryCondition end
struct Clamp <: BoundaryCondition end 

# Neighborhoods
abstract type AbstractNeighborhood{N} end
radius(neighborhood::AbstractNeighborhood) = neighborhood.radius
function cartesians(neighborhood::AbstractNeighborhood{N})::Array{CartesianIndex{N},1} where N
	neighborhood.cartesians
end
(dimensions(::AbstractNeighborhood{N}) where N) = N

struct VonNeumannNeighborhood{N} <: AbstractNeighborhood{N}
	cartesians::Array{CartesianIndex{N},1}
	radius::Int
	function VonNeumannNeighborhood{T}(radius::Int) where T
		arr = CartesianIndex[]
		for i in 1:radius
			for j in 1:T
				zs = zeros(Int, T)
				zs[j] = i
				push!(arr, CartesianIndex(Tuple(zs)))
				zs[j] = -i
				push!(arr, CartesianIndex(Tuple(zs)))
			end
		end
		new(arr, radius)
	end
end

struct MooreNeighborhood{N} <: AbstractNeighborhood{N}
	cartesians::Array{CartesianIndex{N},1}
	radius::Int
	function MooreNeighborhood{T}(radius::Int) where T
		new(
			CartesianIndices(Tuple(-radius:radius for _ in 1:dimensions)) |>
			x -> reshape(x, :),
			radius
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
abstract type AbstractGrid{T,N} <: AbstractArray{T,N} end

function Base.eltype(::Type{AbstractGrid{T,N}}) where {T,N}
	return T
end
(Base.ndims(grid::AbstractGrid{T,N}) where {T,N}) = N
Base.size(grid::AbstractGrid) = size(state(grid))
Base.IndexStyle(::Type{AbstractGrid}) = Base.IndexCartesian

function wrap(::Periodic, grid::AbstractGrid{T,N}, index::CartesianIndex{N})::T where {T,N}
	@inbounds grid[CartesianIndex(mod1.(Tuple(index), size(grid)))]
end
function wrap(::FixedMax, grid::AbstractGrid{T,N}, index::CartesianIndex{N})::T where {T,N}
	maximum(possible_states(grid))
end
function wrap(::FixedMin, grid::AbstractGrid{T,N}, index::CartesianIndex{N})::T where {T,N}
	minimum(possible_states(grid))
end
function wrap(::Clamp, grid::AbstractGrid{T,N}, index::CartesianIndex{N})::T where {T,N}
	@inbounds grid[
		CartesianIndex(
			Tuple(
				clamp(Tuple(index)[i], 1, size(grid)[i]) 
				for i in 1:ndims(grid)
			)
		)
	]
end
@inline function Base.getindex(grid::AbstractGrid{T,N}, index::CartesianIndex{N})::T where {T,N}
    if checkbounds(Bool, grid, index)
       	@inbounds state(grid)[index]
	else
        wrap(boundary_conditions(grid), grid, index)
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

function state(grid::AbstractGrid{T,N})::Array{T,N} where {T,N}
	grid.state
end

function state!(grid::AbstractGrid, newstate)::Nothing
	grid.state = newstate
	return nothing
end
boundary_conditions(grid::AbstractGrid)::BoundaryCondition = grid.bc
(neighborhood(grid::AbstractGrid{T,N})::Array{CartesianIndex{N},1}) where {T,N} = cartesians(grid.neighborhood)
function possible_states(grid::AbstractGrid{T,N})::Array{T,1} where {T,N}
	grid.possible_states
end
function newstate(grid::AbstractGrid{T,N})::Array{T,N} where {T,N}
	similar(state(grid))
end
function neighbors(grid::AbstractGrid{T,N}, i::CartesianIndex{N})::Array{T,1} where {T,N}
    map(neighborhood(grid)) do ind 
		grid[ind + i]
	end
end

function evolve!(grid::AbstractGrid{T,N}; kwargs...)::Nothing where {T,N}
    _newstate = newstate(grid)
    for i in eachindex(grid)
        _newstate[i] = grid(grid[i], neighbors(grid, i); kwargs...)
	end
	state!(grid, _newstate)
	return nothing
end

function tabular_evolution!(grid::AbstractGrid, table::Dict)
    newstate = newstate(grid)
    for i in CartesianIndices(grid)
        newstate[i] = table[grid(neighbors(grid, i))]
	end
	state!(grid, newstate)
	return nothing
end

function simulate!(grid::AbstractGrid{T,N}, max_steps::Int; store_results=true, postrunhook=nothing, kwargs...) where {T,N}
	if store_results 
		results = Array{T,N}[] 
		push!(results, copy(state(grid)))
	end
	if isnothing(postrunhook) 
		postrunhook = (grid, step; kwargs...) -> true
	end
	postrunhook(grid, 0; kwargs...)
	ss = false
	step = 0
	while step < max_steps
		if (!ss) evolve!(grid; kwargs...) end
		if store_results push!(results, copy(state(grid))) end
		message = postrunhook(grid, step; kwargs...)
		ss = message == :shortcurcuit
		if message == :interrupt
			break
		end
		step += 1
	end
	if store_results results else Array{T,N}[] end
end

end
