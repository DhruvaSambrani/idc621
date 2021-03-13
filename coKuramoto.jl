println("Starting script with $(Threads.nthreads())")
using Distributions
using Statistics
using BenchmarkTools
using Plots

println("imported")

#%%
mutable struct Agent
    R::Float64
    V::Float64
end

init(n, σ) = Agent.(1 + 0im, rand(Normal(0., σ), n))
Statistics.mean(arr::Array{Agent}) = mean(cis.(getproperty.(arr, :R)))
function update!(agent::Agent, meanvec, κ, dt)
    agent.R += (agent.V + κ * abs(meanvec) * sin(angle(meanvec) - agent.R)) * dt
end
function evolve!(agents::Array{Agent}, κ::Float64, dt::Float64)
    meanvec = mean(agents)
    Threads.@threads for i in eachindex(agents)
        update!(agents[i], meanvec, κ, dt)
    end
    meanvec
end
function simulate(steps::Int, n::Int, σ::Float64, κ::Float64, dt::Float64)
    agents = init(n, σ)
    means = Array{Complex,1}(undef, Int(steps))
    for step in 1:steps
        means[step] = evolve!(agents, κ, dt)
    end
    means
end


#%%

# @show @benchmark evolve(init(500, 0.3), 0.5, 1e-3)

#%%


function graph(κs, σs)
    arr = zeros(Float64, length(κs), length(σs))
    for (i, κ) in enumerate(κs), (j, σ) in enumerate(σs)
        println(κ, " ", σ)
        arr[i,j] = mean(abs, simulate(Int(1e5), 500, σ, κ, 1e-4)[end - 1000:end])
    end
    arr
end
# %%

κs, σs = 0:0.1:2, 0:0.1:2
zs = graph(κs, σs)
# %%
plot(σs, κs, zs, st=:surface, zlims=(0, 1), xlabel="Sigma", ylabel="Kappa", zlabel="r_inf", title="R∞(σ, κ)")
savefig("ksr")

# %%
# plot(abs.(simulate(Int(1e5), 500, 1., 0., 1e-4)))
