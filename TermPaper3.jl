### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 8d7c0eb0-84ca-11eb-0018-bdb2e389ee52
begin
	using PlutoUI
	using Distributions
	using Statistics
	using Plots
end

# ╔═╡ a34a7c40-84ca-11eb-3433-41615dfc1c08
mutable struct Agent
	R::Float64
	V::Float64
end

# ╔═╡ acefa9a0-84ca-11eb-1f9d-4918cbdd1211
init(n, σ) = Agent.(1., rand(Normal(0., σ), n))

# ╔═╡ c12bde70-84ca-11eb-1fe2-634e16d69909
Statistics.mean(arr::Array{Agent})::Complex = mean(map(arr) do ag; cis(ag.R); end)

# ╔═╡ c4b72810-84ca-11eb-1855-23703a18b3e4
function update!(agent::Agent, κr::Float64, Phi::Float64, dt::Float64)
	agent.R += (agent.V + κr * sin(Phi - agent.R)) * dt
end

# ╔═╡ cb3e73a0-84ca-11eb-3c68-8fc3fc5cf239
function evolve!(agents::Array{Agent}, κ::Float64, dt::Float64)
	meanvec = mean(agents) 
	update!.(agents, κ*abs(meanvec), angle(meanvec), dt)
	meanvec
end

# ╔═╡ d00b64fe-84ca-11eb-1e68-f73bbca9e373
function simulate(steps::Int, n::Int, σ::Float64, κ::Float64, dt::Float64)
	agents = init(n, σ)
	means = Array{Complex,1}(undef, Int(steps))
	for step in 1:steps
		means[step] = evolve!(agents, κ, dt)
	end
	means
end

# ╔═╡ 8c9b5c8e-84d3-11eb-0f9d-47c844a46a68
means = simulate(100000, 1000, 1., 2., 5e-4);

# ╔═╡ c618ba40-84d2-11eb-3990-93a786c8f8a9
plot(abs.(means))

# ╔═╡ Cell order:
# ╠═8d7c0eb0-84ca-11eb-0018-bdb2e389ee52
# ╠═a34a7c40-84ca-11eb-3433-41615dfc1c08
# ╠═acefa9a0-84ca-11eb-1f9d-4918cbdd1211
# ╠═c12bde70-84ca-11eb-1fe2-634e16d69909
# ╠═c4b72810-84ca-11eb-1855-23703a18b3e4
# ╠═cb3e73a0-84ca-11eb-3c68-8fc3fc5cf239
# ╠═d00b64fe-84ca-11eb-1e68-f73bbca9e373
# ╠═8c9b5c8e-84d3-11eb-0f9d-47c844a46a68
# ╠═c618ba40-84d2-11eb-3990-93a786c8f8a9
