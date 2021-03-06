using Distributions
using Statistics
using JLD

println("imported")
# %%
mutable struct Agent
    R::Float64
    V::Float64
end

init(n, σ) = Agent.(1., rand(Normal(0., σ), n))
Statistics.mean(arr::Array{Agent}) = mean(map(arr) do a; cis(a.R); end)
function update!(agent::Agent, r::Float64, Ψ::Float64, κ::Float64, dt::Float64)
    agent.R += (agent.V + κ * r * sin(Ψ - agent.R)) * dt
end
function evolve!(agents::Array{Agent}, κ::Float64, dt::Float64)
    meanvec = mean(agents)
    r, Ψ = abs(meanvec), angle(meanvec)
    update!.(agents, r, Ψ, κ, dt)
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

function graph(κs::Array{Float64}, σs::Array{Float64})
    arr = zeros(Float64, length(κs), length(σs))
    for (i, κ) in enumerate(κs), (j, σ) in enumerate(σs)
        arr[i,j] = mean(
                abs, 
                @view simulate(Int(1e5), 1000, σ, κ, 1e-4)[end - 1000:end]
	)
    end
    arr
end
# %%
κs, σs = collect(0.0:0.1:2.5), collect(0.0:0.1:2.5)
zs = graph(κs, σs)
# %%
save("save.jld", "zs", zs, "ks", κs, "ss", σs)

