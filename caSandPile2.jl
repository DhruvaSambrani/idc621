@time using Plots
@time using LsqFit
@time using Statistics
# %%
@inline function getindexbc(grid, i)
    if checkbounds(Bool, grid, i)
        @inbounds return grid[i]
    else
        return 0
    end
end
@inline iscritical(i::Int, hc::Int) = i >= hc
getindexbc(grid) = i -> getindexbc(grid)
function evolveat(
        grid::Array{Int,2}, 
        hood::Array{CartesianIndex{2},1}, 
        index::CartesianIndex, 
        hc::Int,
        topple::Array{Int},
        step::Int
    )
    res = grid[index]
    if res >= hc
        res -= 4
        @inbounds topple[step] += 1
    end
    c = count(h -> iscritical(getindexbc(grid, h + index), hc), hood)
    res = res + c
    return res
end
function evolve!(
        grid::Array{Int,2}, 
        simgrid::Array{Int,2},
        hood::Array{CartesianIndex{2},1},
        hc::Int,
        topple::Array{Int,1},
        step::Int
    )
    for index in CartesianIndices(grid)
        @inbounds simgrid[index] = evolveat(grid, hood, index, hc, topple, step)
    end
end
function ready_next!(grid::Array{Int,2}, hc::Int, topple::Array{Int,1}, step::Int, max_steps::Int)
    if count(iscritical.(grid, hc)) == 0
        if step == max_steps
            return false, step
        else
            xrand, yrand = rand.(Base.OneTo.(size(grid)))
            grid[xrand, yrand] += 1
            step += 1
        end
    end
    return true, step
end

function simulate(grid::Array{Int}, hc::Int, max_steps::Int)
    hood = CartesianIndex.([(-1, 0), (1, 0), (0, -1), (0, 1)])
    topple = fill(0, max_steps)
    step = 1
    st = true
    simgrid = similar(grid)
    results = Array{Int}[]
    while st
        evolve!(grid, simgrid, hood, hc, topple, step)
        grid, simgrid = simgrid, grid
        st, step = ready_next!(grid, hc, topple, step, max_steps)
    end
    return results, topple
end
# %%

hc = 4
results, topple = simulate(rand(1:hc - 1, 20, 20), 4, 10000);

# %%

function freqtable(k::Array{T}) where T
    d = Dict{T,Int}()
    for i ∈ k
        if i != 0
            d[i] = get(d, i, 0) + 1
        end
    end
    return d
end
function getlogdata(topple_freq)
    x = log.(keys(topple_freq))
    perm = sortperm(x)
    y = log.(values(topple_freq))[perm]
    x = x[perm]
    return x, y
end
function bin(xs, ys, length)
    binx = range(x[1], stop=x[end], length=length)
    _step = step(binx)
    biny = [
        mean(ys[findall(x -> val <= x < (val + _step), xs)])
        for val in binx
    ]
    binx = (binx .+ _step)[1:end - 1]
    return binx, biny
end

x, y = getlogdata(freqtable(topple))
binx, biny = bin(x, y, 100)
plot(binx, biny, st=:scatter)
# %%
@. model(x, p) = p[1] * x + p[2]
fit = curve_fit(model, x[1:45], y[1:45], [0.,0.])
println(fit.param)
xp = 0:0.01:maximum(x)
yp = model(xp, fit.param)
plot!(xp, yp)
