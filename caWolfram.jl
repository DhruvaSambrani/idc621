#%%
include("CellularBase.jl")
using .CellularBase
import .CellularBase: evolve!
using Plots
using Images
import StatsBase: sample
#%%
mutable struct WolframGrid <: AbstractGrid
	state::BitArray
	bc::BoundaryCondition
	neighborhood::AbstractNeighborhood
	possible_states::Array
	function WolframGrid(arr::BitArray)
		new(arr, Periodic, MooreNeighborhood(1, 1), [true, false])
	end
	function WolframGrid(size::Int)
		s = falses(size)
		s[size ÷ 2] = true
		WolframGrid(s)
	end
	function WolframGrid(trueratio::Float64, size::Int)
		ts = sample(1:size, ceil(Int, trueratio * size), replace=false)
		s = falses(size)
		s[ts] .= true
		WolframGrid(s)
	end
end

function (g::WolframGrid)(neighbors)
	sum = 0
	for (i, v) in enumerate(reverse(neighbors))
		sum += v * (2^(i - 1))
	end
	sum
end

function wr_table(rule::Int)
	Dict{Int,Int}([(i, rule >> i & 1) for i in 0:7])
end

function CellularBase.evolve!(grid::WolframGrid; kwargs...)	
	tabular_evolution!(grid, wr_table(kwargs[:rule]))
end

function makeimg(k)
	k .|> (x -> x .== false) |> x->hcat(x...)' .|> Gray
end
#%%
k = simulate!(WolframGrid(100), 100, rule=26);
plot(makeimg(k))
#%%
