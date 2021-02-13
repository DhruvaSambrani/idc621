#%%
push!(LOAD_PATH, ".");
#%%
using CellularBase
using Plots
#%%
mutable struct SandPile<:AbstractGrid{Int}
    state::Array{Int, 2}
    otherstate::Array{Int, 2}
    h_c::Int
    neighborhood::AbstractNeighborhood
    bc::BoundaryCondition
    function SandPile(state::Array{Int, 2}, h_c::Int)
        new(state, copy(state), h_c, VonNeumannNeighborhood(1,2), FixedMin)
    end
end

CellularBase.possible_states(sp::SandPile) = 0:sp.h_c

CellularBase.newstate(grid::SandPile) = grid.otherstate;
function CellularBase.state!(grid::SandPile, newstate)
    grid.otherstate = grid.state
    grid.state = newstate
end

iscritical(val::Int, g::SandPile) = val >= g.h_c
iscritical(g::SandPile) = p -> iscritical(p, g)

function (_g::SandPile)(neighbors::Array{Int}; kwargs...)
    current = neighbors[1]
    retval = current
    if iscritical(current, _g)
        retval -= _g.h_c
        kwargs[:topple][end] += 1
    end
    retval += count(iscritical(_g), neighbors[2:end])
    return retval
end

function ready_next(grid::SandPile, step::Int; topple::Array{Int}, max_steps)
    if count(iscritical(grid), state(grid)) == 0
        # No critical left
        if length(topple) == max_steps
            # Enough, my god
            return :interrupt
        else
            # Moar, moar, MOAR!
            # xrand, yrand = rand.(Base.OneTo.(size(grid)))
            xrand, yrand = 5,5
            # Potentially destabilize - HAHAHAA
            grid.state[xrand, yrand] += 1
            # Make new counter
            push!(topple, 0)
        end
    end
end
#%%
h_c = 4
topple_count = Int[]
linear(i,j) = trunc(Int, (h_c-1)*(1 - √(i^2+j^2)/√50))

arr = [linear(5-i,5-j) for i in 1:9, j in 1:9]

grid = SandPile(zeros(Int, 9,9), h_c)
simulate!(grid, typemax(Int), postrunhook=ready_next, topple=topple_count, max_steps=1000)

topple_freq = freqtable(topple_count)

plot(keys(topple_freq), values(topple_freq))
#=
#%%
anim = @animate for frame ∈ results
    wireframe(frame, zlims=(0,4))
end
gif(anim, fps=30)
#%%=#

function freqtable(k::Array{T}) where T
    d = Dict{T, Int}()
    for i ∈ k
        d[i] = get(d, i, 0) + 1
    end
    return d
end