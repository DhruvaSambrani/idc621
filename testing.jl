# %%
if !("." in LOAD_PATH) push!(LOAD_PATH, ".") end
using Revise
@time using CellularBase
using BenchmarkTools
using Random
# %%
mutable struct SandPile <: AbstractGrid{Int,2}
    state::Array{Int,2}
    otherstate::Array{Int,2}
    h_c::Int
    neighborhood::AbstractNeighborhood
    bc::BoundaryCondition
    function SandPile(state::Array{Int,2}, h_c::Int)
        new(state, similar(state), h_c, VonNeumannNeighborhood{2}(1), FixedMin())
    end
end
@inline iscritical(c::Int, g::SandPile) = c >= g.h_c
@inline iscritical(g::SandPile) = i -> iscritical(i, g)

CellularBase.newstate(s::SandPile) = s.otherstate
function CellularBase.state!(s::SandPile, newstate) 
    temp = s.state
    s.state = s.otherstate
    s.otherstate = temp
end

function (_g::SandPile)(current::Int, neighbors::Array{Int, 1})
    h = iscritical(_g)
    return current + h(current)*-4 + count(h, neighbors)
end

CellularBase.possible_states(g::SandPile) = 0:(g.h_c + 3)
# %%

function getindexbc(grid, i)
    if checkbounds(Bool, grid, i)
        @inbounds return grid[i]
    else
        return 0
    end
end
getindexbc(grid) = i -> getindexbc(grid)
function evolveat(grid, hood, index, hc)
    res = grid[index]
    if res >= hc
        res -= 4
    end
    c = count(h -> hc <= getindexbc(grid, h + index), hood)
    res = res + c
    return res
end
function evolve(grid::Array{Int}, hood::Array{CartesianIndex{2},1}, hc::Int; kwargs...)
    simgrid = similar(grid)
    for index in CartesianIndices(grid)
        simgrid[index] = evolveat(grid, hood, index, hc)
    end
    return simgrid
end

function baresimulate(grid::Array{Int}, hc::Int, steps::Int)
    hood = CartesianIndex.([(-1, 0), (1, 0), (0, -1), (0, 1)])
    for i in 1:steps
        grid = evolve(grid, hood, hc)
    end
    return grid
end
function baresimulate!(cbgrid::SandPile, steps)
    for i in 1:steps
        evolve!(cbgrid)
    end
end

function gridtest(cbgrid::SandPile)
    for i = 1:400
        cbgrid(1, [1, 1, 1, 1])
    end
end
function gridtest(ngrid::Array{Int})
    hood = cartesians(VonNeumannNeighborhood{2}(1));
    for i in 1:400
        evolveat(ngrid, hood, CartesianIndex(1,1), 4)
    end
end

CellularBase.boundary_conditions(p::SandPile)::FixedMin = FixedMin()

# %%
println("===============")
Random.seed!(0);
hc = 4;
ngrid = rand(1:hc, 20, 20);
cbgrid = SandPile(copy(ngrid), 4);
hood = cartesians(VonNeumannNeighborhood{2}(1));

@btime baresimulate(ngrid, 4, 5)
# @btime evolve(ngrid, hood, 4)
# @btime evolve(cbgrid)
# @btime CellularBase.evolve!(cbgrid)
@btime simulate!(cbgrid, 5, store_results=false)
# @btime baresimulate!(cbgrid, 1)

# %%
