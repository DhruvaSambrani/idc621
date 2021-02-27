# %%
push!(LOAD_PATH, ".");
# %%
@time using CellularBase
#@time using Plots
#@time using LsqFit
@time using Random
# %%
mutable struct SandPile <: AbstractGrid{Int8, 2}
    state::Array{Int8,2}
    h_c::Int8
    neighborhood::AbstractNeighborhood
    bc::BoundaryCondition
    function SandPile(state::Array{Int8,2}, h_c::Int8)
        new(state, h_c, VonNeumannNeighborhood(1, 2), FixedMin)
    end
end

CellularBase.possible_states(sp::SandPile) = 0:sp.h_c

@inline iscritical(val::Int8, g::SandPile) = val >= g.h_c
@inline iscritical(g::SandPile) = p -> iscritical(p, g)

function (_g::SandPile)(current::Int8, neighbors::Array{Int8}; kwargs...)
    h = iscritical(_g)
    if h(current) kwargs[:topple][end] += 1 end
    return current + (h(current) ? -4 : 0) + count(h, neighbors)
end

function ready_next(grid::SandPile, step::Int; topple::Array{Int64}, max_steps)
    if count(iscritical(grid), state(grid)) == 0
        # No critical left
        if length(topple) == max_steps
            # Enough, my god
            return :interrupt
        else
            # Moar, moar, MOAR!
            xrand, yrand = rand.(Base.OneTo.(size(grid)))
            # Potentially destabilize - HAHAHAA
            grid.state[xrand, yrand] += 1
            # Make new counter
            push!(topple, 0)
        end
    end
    return :continue
end
# %%
Random.seed!(0)
h_c = Int8(4)
topple_count = Int[0]
w = 20
rand_arr = rand(Int8.(1:h_c-1), w, w)
zero_arr = zeros(Int8, w, w)

grid = SandPile(rand_arr, h_c)
@time simulate!(
    grid, typemax(Int), 
    postrunhook=ready_next, topple=topple_count,
    max_steps=100_000, store_results=false);

# %%

function freqtable(k::Array{T}) where T
    d = Dict{T,Int}()
    for i ∈ k
        d[i] = get(d, i, 0) + 1
    end
    return d
end


topple_freq = freqtable(topple_count)

topple_freq = filter!((p) -> p[1] * p[2] != 0, topple_freq)

x = log.(keys(topple_freq))
y = log.(values(topple_freq))

plot(x,y, st=:scatter)
# %%
@. model(x, p) = p[1] * x + p[2]
fit = curve_fit(model, x, y, [0.,0.])
println(fit.param)
xp = 0:0.01:maximum(x)
yp = model(xp, fit.param)
plot!(xp, yp)

# %%
#=
anim = @animate for frame ∈ results
    wireframe(frame, zlims=(0,4))
end
gif(anim, fps=30)
=#
#%%
