using Distributions
using Statistics
using Plots

println("imported")
#%%
struct Agent
    R::Float64
    V::Float64
end

init(n, σ) = Agent.(1 + 0im, rand(Normal(0., σ), n))
Statistics.mean(arr::Array{Agent}) = mean(cis.(getproperty.(arr, :R)))
function update(agent::Agent, meanvec, κ, dt) 
    Agent(
        agent.R + (
                agent.V + 
                κ * abs(meanvec) * sin(angle(meanvec) - agent.R)
            ) * dt,
        agent.V
    )
end
function evolve(agents::Array{Agent}, κ, dt)
    meanvec = mean(agents)
    meanvec, update.(agents, meanvec, κ, dt)
end
function simulate(steps::Int, n::Int, σ::Float64, κ::Float64, dt::Float64)
    agents = init(n, σ)
    means = Array{Complex,1}(undef, Int(steps))
    records = Array{Complex,1}(undef, Int(steps))
    for step in 1:steps
        means[step], agents = evolve(agents, κ, dt)
    end
    means
end
function graph(κs, σs)
    arr = zeros(Float64, length(κs), length(σs))
    for (i, κ) in enumerate(κs), (j, σ) in enumerate(σs)
        println(κ, " ", σ)
        arr[i,j] = mean(abs, simulate(Int(1e5), 500, σ, κ, 1e-4)[end - 1000:end])
    end
    arr
end
# %%
κs, σs = 0:0.2:2, 0:0.1:1
zs = graph(κs, σs)
# %%
plotlyjs()
plot(σs, κs, zs, st=:surface, zlims=(0, 1), xlabel="Sigma", ylabel="Kappa", zlabel="r_inf")
# %%
plot(abs.(simulate(Int(1e5), 500, 1., 0., 1e-4)))
