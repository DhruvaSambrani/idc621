### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 1bc961f0-66f7-11eb-17c8-2b54d4c4328d
begin
	push!(LOAD_PATH, ".")
	using CellularBase
end

# ╔═╡ f8149ae0-66f6-11eb-0335-39baafdb24c8
begin
	using Plots
	using Images
	using PlutoUI
	using Statistics
end

# ╔═╡ 9a74ef90-67ad-11eb-1fc3-e7080fbfc7e1
html"""<center><h1 style="font-size: 3.5rem">Term Paper 1<h1></center><br>
This work is a submission for the course Modelling Complex Systems(IDC-631) for the year 2020-21 Spring (Even) Semester by <b>Dhruva Sambrani, MS18163</b>
"""

# ╔═╡ eb2aad32-687f-11eb-0c6b-91b443cdb3d3
md"""# Abstract
This paper looks at the effect of a Lockdown on the spread of a pathogen. This setup is modelled via a Probabilisitic/Non Deterministic 2D cellular automaton.

This paper first [introduces](#7cb22830-67c9-11eb-0e0e-41f2316ed468) some working definitions along with the same definitions in code.
Then features are added, and parallelly introduced to the reader, [layer](#c0318a00-67ac-11eb-39b7-e3c13a5aff57) upon [layer](#1cf29850-66fa-11eb-3f76-db755c0a4c50) with logical reasoning for their implementations.
Following this, the actual topic of the paper is introduced with some [results](#5b1e6e50-6877-11eb-1ca4-3df17b85575e) that show the efficiency of [Lockdowns](#9e71a7c0-67f1-11eb-0282-672d358d2292).
Finally, [future scope of this work](#c7ea4fa0-687b-11eb-0dbc-7f316bd4aefd) is also discussed.

!!! note "Note to Evaluator"
    The author believes that papers which involve a lot of simulation and code should also be accompanied with the code they use so that readers can gain deeper insight into the processes described in the paper. Hence, the author has attempted to create a paper which interlaces the information with the corresponding code.

The coding language used is `Julia`[^julia] and following document is made with `Pluto.jl` [^pluto].

"""

# ╔═╡ 6f7bce82-67ad-11eb-1ec1-79ce2620fe9b
PlutoUI.TableOfContents(title="Term Paper 1 - Index")

# ╔═╡ 3459c040-6858-11eb-09c1-95c7bdb02722
html"<style>
.plutoui-toc.aside {
    top: 15%;
    bottom: 20%;
}</style>"

# ╔═╡ 7cb22830-67c9-11eb-0e0e-41f2316ed468
md"# Code Setup"

# ╔═╡ 65ff1d80-66f9-11eb-11c2-69c47b8458e4
md"## `Include packages`"

# ╔═╡ 330799c0-6884-11eb-230d-71c512e4a447
md"`CellularBase`[^cellularbase] is a package written by the author themselves to abstract out certain repetitive funcitions for all Cellular Automata."

# ╔═╡ db4fd9e0-6883-11eb-384e-6df2c1d7a52b
md"The `Plots.jl`[^plots] and `Images.jl`[^images] packages are used for plotting and the `Statistics`[^stats] Package for some statistical analysis. `PlutoUI.jl`[^plutoui] is used for some UI elements such as the ToC"

# ╔═╡ a95c8c70-67ad-11eb-124b-7b66d5547321
md"## Define `Structures` and `Functions`"

# ╔═╡ 73ea0ae0-66f9-11eb-19cf-b7ab6151ef63
md"### Define `LifeGrid`"

# ╔═╡ 038cbdf2-6880-11eb-1b2b-2b4c4974be49
md"""This is the basic structure that defines all the parameters of a particular "grid". The name `LifeGrid` was given by mistake but was not changed as it would entail too many other changes"""

# ╔═╡ 35d0be40-66f7-11eb-2f60-01982265ccb0
mutable struct LifeGrid <: AbstractGrid{Int}
	state::Array{Int}
	otherstate::Array{Int}
	bc::BoundaryCondition
	neighborhood::AbstractNeighborhood
	τi::Int64
	τr::Int64
	max::Int64
	function LifeGrid(
			size::Tuple{Int,Int}, τi::Int, τr::Int;
			max::Int64, init=nothing, 
			nh::AbstractNeighborhood=VonNeumannNeighborhood(1, 2))
		z = zeros(size)
		if isnothing(init) z[(size .÷ 2)...] = 1
		else z[rand(CartesianIndices(z), init)] .= 1
		end
		new(z, copy(z), FixedMax, nh, τi, τr, max)
	end
	function LifeGrid(
			arr::Array{Int}, τi::Int, τr::Int;
			max::Int64,
			nh::AbstractNeighborhood=VonNeumannNeighborhood(1, 2)
		)
		new(arr, copy(arr), FixedMax, nh, τi, τr, max)
	end
end

# ╔═╡ 91afcec0-66f9-11eb-37b1-f78de5e0ba8a
md"### Define method of evolution"

# ╔═╡ fcc23790-686e-11eb-03d7-1f0946446f49
md"The method of evolution follows the following algorithm -
1. If the agent is susceptible, infect the agent depending on neighborhood and return. The exact dependacy is discussed later in [Non Deterministic runs](#2563e210-67d5-11eb-089e-f376da973765)
2. Otherwise, the state is increased by 1
3. If the state is recovered, set it to 0

Mathematically this can be written as

$s_{i+1} = \begin{cases}
    1, & s_i=0 \text{ and } \rho(\text{neighbours}) \\
    0, & s_i=0 \text{ and } 1-\rho(\text{neighbours}) \\
	(s_i+1) \mod (\tau_i+\tau_r+1) & \text{otherwise}
\end{cases}$

"

# ╔═╡ b8e7c920-66f9-11eb-397f-3760662a648c
md"### Define States"

# ╔═╡ 8b6833e0-6870-11eb-2593-69edaec0a7b8
md"States are defined in the following way
1. Susceptible $\implies 0$
2. Infected $\implies i \in [1, \tau_i]$
3. Refactory $\implies i \in [\tau_i + 1, \tau_i+\tau_r]$
4. Recovered $\implies i = (\tau_i+\tau_r + 1)$
"

# ╔═╡ 21d1eb20-66f8-11eb-26d3-4fb244ca8646
begin
	infect(_count::Int64, g::LifeGrid) = if (_count == 0) false else rand() < 1 / g.max * _count end
	susceptible(i::Int, g::LifeGrid) = i == 0
	infectious(i::Int, g::LifeGrid) = i ∈ 1:g.τi
	refactory(i::Int, g::LifeGrid) = i ∈ (g.τi + 1):(g.τr + g.τi)
	recovered(i::Int, g::LifeGrid) = i == (g.τr + g.τi) + 1
end

# ╔═╡ 263db440-66f9-11eb-03b5-1fce784fad2d
function (_g::LifeGrid)(neighbors::Array{Int}; kwargs...)
	current = neighbors[1]
	if susceptible(current, _g)
		if infect(count(i -> infectious(i, _g), neighbors), _g)
			return 1
		else
			return 0
		end
	else
		return (current + 1) % (_g.τi + _g.τr)
	end
end

# ╔═╡ de865ed2-66f9-11eb-31dc-c782d4338af3
md"### Define helper `functions`"

# ╔═╡ c3e6cc10-66f7-11eb-039f-df82addf1a5e
begin
	CellularBase.possible_states(l::LifeGrid) = 0:(l.τi + l.τr)
	CellularBase.newstate(grid::LifeGrid) = grid.otherstate;
	function CellularBase.state!(grid::LifeGrid, newstate)
	    grid.otherstate = grid.state
	    grid.state = newstate
	end
end

# ╔═╡ 4d4cdcc0-67b0-11eb-2278-8de64128d6f1
function Base.show(io::IO, ::MIME"text/html", grid::LifeGrid)
	final = Markdown.parse(
		"---" *
		"\n\n**Grid Type:** " * repr(typeof(grid), context=io) *
		"\n\n**\$\\tau_i:\$** " * string(grid.τi) *
		"\n\n**\$\\tau_r:\$** " * string(grid.τr) *
		"\n\n**Max neighbors for probabilistic spread:** " * string(grid.max) *
		"\n\n**Possible States:** " * string(possible_states(grid)) *
		"\n\n**Boundary Condition:** " * string(boundary_conditions(grid)) *
		"\n\n**Size:** " * string(size(grid)) *
		"\n\n**Neighborhood:** "
	)
	show(io, MIME("text/html"), final)
	show(io, MIME("text/html"), grid.neighborhood)
	show(io, MIME("text/html"), md"---")
end

# ╔═╡ 7cdebc60-684b-11eb-38a2-45873e4b51d0
function makegraph(ldgrid, ldresults, ldradii)
    p = plot(1:length(ldresults), ldresults .|> frame -> count(frame .|> p -> susceptible(p, ldgrid)) / length(state(ldgrid)), color=:green, label="Sus")
    plot!(1:length(ldresults), ldresults .|> frame -> count(frame .|> p -> infectious(p, ldgrid)) / length(state(ldgrid)), color=:red, label="Inf")
    plot!(1:length(ldresults), ldresults .|> frame -> count(frame .|> p -> refactory(p, ldgrid)) / length(state(ldgrid)), color=:blue, label="Ref")
    q = plot(ldradii, label="Neighbourhood Radius", legend=:bottomright)
    plot(p, q, size=(600, 500), layout=(2, 1))
end;

# ╔═╡ 5031733e-66f9-11eb-107f-6d69089d9384
function color(i, grid)
    if susceptible(i, grid)
        RGB(0.84, 0.95, 0.80)
    elseif infectious(i, grid)
        RGB(0.83, 0.11, 0.31)
    else
        RGB(0.16, 0.11, 0.42)
    end
end;

# ╔═╡ c0318a00-67ac-11eb-39b7-e3c13a5aff57
md"# Deterministic Case - A First Run"

# ╔═╡ 47b81340-67cb-11eb-372e-972e4d4ff746
md"## Introduction
As a first trial, we can try to run the simulation for the deterministic case, where the infection can spread even if 1 neighbor is infected. We also verify that the results observed are very similar to results obtained in [^Moitra2019]
"

# ╔═╡ 02c24a70-67ae-11eb-2944-61a7963eb947
LifeGrid((100, 100), 4, 9, max=1, nh=VonNeumannNeighborhood(1, 2))

# ╔═╡ fbdba862-67cf-11eb-39bc-f37eb3662785
md"## `state(center patient) = 1`"

# ╔═╡ fe0dda70-66fa-11eb-0e02-73c05666aace
dgrid, dresults = let
	grid = LifeGrid(
		(100, 100), 
		4, 9,
		max=1,
		nh=VonNeumannNeighborhood(1, 2));
	results = simulate!(grid, 150)
	grid, results
end;

# ╔═╡ 8315d810-677a-11eb-3953-8faad42ada40
@gif for frame ∈ dresults[5:end]
	plot(frame .|> p -> color(p, dgrid))
end

# ╔═╡ 7d24e300-67d0-11eb-0864-41c6166f13ac
md"## Random states from 0 to $\tau_i$"

# ╔═╡ b8523010-67ce-11eb-2396-4780284dfa01
randgrid, randresults = let
	arr = rand(0:4, 100, 100)
	grid = LifeGrid(
		arr,
		4, 5,
		max=1,
		nh=VonNeumannNeighborhood(1, 2));
	results = simulate!(grid, 50)
	grid, results
end;

# ╔═╡ cf22a450-67ce-11eb-2141-fd6f3efc689a
@gif for frame ∈ randresults
	plot(frame .|> p -> color(p, randgrid))
end

# ╔═╡ 07c440a0-67d1-11eb-3a8b-d559cac6276e
md"## `state(center agent) = τᵢ; τᵢ=τᵣ`"

# ╔═╡ 30ce7a10-67bd-11eb-027d-8dc583201bc1
ϕgrid, ϕresults = let
	arr = zeros(Int, 100, 100)
	arr[50,50] = 4
	grid = LifeGrid(
		arr,
		4, 4,
		max=1,
		nh=VonNeumannNeighborhood(1, 2));
	results = simulate!(grid, 50)
	grid, results
end;

# ╔═╡ 4ca7fe50-67bd-11eb-382c-1304f3c1c3c4
@gif for frame ∈ ϕresults
	plot(frame .|> p -> color(p, ϕgrid))
end

# ╔═╡ 643ffb90-67ad-11eb-20c6-7744050e5448
md"""## Observations in the deterministic case

We can see that the infection spreads outward and then dies off.

The refactory agents act as a buffer between the leading edge of infected agents and the recovered(susceptible) agents inside the wave. This is similar to how a batch of burnt material doesn't allow the fire to move back in.

!!! tip "Spread Back"
    The infection can spread back into the wave only if either - 
    -
      1. the radius of the neighborhood is larger than the refactory period 
      2. if $\tau_i \ge \tau_r$ and $\text{state(patient zero)}\ge\tau_r$


"""

# ╔═╡ 1cf29850-66fa-11eb-3f76-db755c0a4c50
md"# Non Deterministic Case"

# ╔═╡ 393a5bc0-67cb-11eb-3be1-a9ec7943b5f3
md"""## Why?
Mainly because determinism is so 1900s.

Seriously though, the reason why we don't see too many interesting things in the deterministic case is because of the refactory buffer separating the infected leading edge and the recovered members.

Allowing some agents to escape the infection wave allows breaks in the refactory buffer allowing the infection to come back into the recovered agents. It seems very sadistic to play for the pathogen, but what's the fun otherwise.
"""

# ╔═╡ 2563e210-67d5-11eb-089e-f376da973765
md"""## How is it implemented?
In code, obviously. Oh you mean what's the logic behind it? 

If the spread is probabilistic, one way to do it would be to set a non zero probability if the number of infected neighbors is non-negative. While this does make the spread probabilistic, it doesn't model _real life_ very well. In real life, one's chances of infection increases with the increase in infected neighbors. While the exact distribution is unknown, we can model the spread by assuming that the chance of spreading follows 

$\rho_\text{infect}{}(\text{count}) = \begin{cases}
  \frac{\text{count}}{\text{max}} & 0\leq \text{count}\leq \text{max} \\
  1 & \text{max}\leq \text{count} 
\end{cases}$

where $\text{max}$ is a parameter of the simulation. This basically means that for every infected neighbour, the chance of infection is linear, until a limit, after which the agent will always get infected.

This definition has the added benefit of being independent on the definition of the neighbourhood, and can be attributed only to the pathogen in question. 

Every `LifeGrid` has a critical limit of infected neighbors, denoted by `max`. For every agent, $\rho$ is calculated and infection is done by a Bernoulli trial.

!!! info "High values of max"

    Hence, it also makes sense to set `max` higher than the size of the neighbourhood itself, for that would allow us to set the highest possible probability to be lesser than 1. 
    
    After all, what if the neighborhood changes during the run? How? Keep reading
"""

# ╔═╡ dafac300-67e9-11eb-0972-f9070193399e
md"## Max = 4"

# ╔═╡ 46055e90-67e9-11eb-2b6f-dffa58138005
ρgrid, ρresults = let
	grid = LifeGrid(
		(100, 100),
		4, 5,
		max=4,
		nh=VonNeumannNeighborhood(1, 2));
	results = simulate!(grid, 400)
	grid, results
end;

# ╔═╡ 5626db00-67e9-11eb-0e12-a77b6b9c58be
@gif for frame ∈ ρresults
	plot(frame .|> p -> color(p, ρgrid))
end every 2

# ╔═╡ d4c2cff0-67e9-11eb-1f5b-03c7bab91782
md"## Max = 9"

# ╔═╡ 3056cdae-67ec-11eb-0b35-b3d6cde9b0fa
ρ7grid, ρ7results = let
	grid = LifeGrid(
		(100, 100),
		4, 5,
		max=9,
		init=5,
		nh=VonNeumannNeighborhood(1, 2));
	results = simulate!(grid, 100)
	grid, results
end;

# ╔═╡ 3b0024a0-67ec-11eb-3398-eb07250758a9
@gif for frame ∈ ρ7results
	plot(frame .|> p -> color(p, ρ7grid))
end

# ╔═╡ 6d822d52-67ed-11eb-1b71-679d491bf133
md"## Max = 2"

# ╔═╡ 76e43730-67ed-11eb-201f-53bf8d3e3bfb
ρ2grid, ρ2results = let
	grid = LifeGrid(
		(100, 100),
		4, 5,
		max=2,
		nh=VonNeumannNeighborhood(1, 2));
	results = simulate!(grid, 150)
	grid, results
end;

# ╔═╡ 8445d050-67ed-11eb-2411-f776a43c7f40
@gif for frame ∈ ρ2results
	plot(frame .|> p -> color(p, ρ2grid))
end every 2

# ╔═╡ 0f5ba100-67ef-11eb-15f5-3fbaad70dde2
md"""## Observations

Evidently, higher max values imply smaller probability of spreading. This is evident in the run with `max = 9`, where very few infected cases are seen at any point in time, and agents often recover before they spread it on, leading to the pathogen dying out. 

On the flip side, if `max` is too close to 1, then there it acts very similar to the deterministic case _most of the time_. This is also evident in the fact that the case where `max` is 1, it is exactly the deterministic case which we've seen also dies out.

!!! tip "Relevant Quote"
    It is possible to have too much of a good thing
    
    \- **Aesop**
"""

# ╔═╡ 9e71a7c0-67f1-11eb-0282-672d358d2292
md"""# Lockdown
## Introduction
The word _lockdown_ seems to have taken over our lives and completely changed the way that we go about our daily lives. But why are we doing this? And how effective are these lockdowns? We explore this question by looking at how we can simulate a lockdown period for our agents and see how this affects the infection spread.
"""

# ╔═╡ ef666b40-67f4-11eb-2f31-05145a6dfcba
md"""## The logic behind lockdowns

!!! info "Definition"
    A lockdown is a restriction policy for people or community to stay where they are, usually due to specific risks to themselves or to others if they can move and interact freely - [Wikipedia](https://en.wikipedia.org/wiki/Lockdown)

From the definition, it is clear that the intent is to prevent people from interacting and hence hopefully reducing the spread of the disease. If an agent interacts with fewer people, they have a lower chance of getting infected.
"""


# ╔═╡ 5586a260-67f4-11eb-11e2-b706fd0d4dfe
md"""## Neighborhoods, Virulence and Infection Rate
Neighbourhood size and intrinsic virulence of a pathogen together decide the pathogen's spreading capabilities. 

Neighbourhood size is a representation of the amount of inter-agent interactions. Due to the way the model is defined, it also shows how _far_ the agent can influence. It defines the number of agents one agent can meet.

Intrinsic virulence of a pathogen is how well it spreads given that two agents meet. There is sufficient discussion regarding this in the above section.

If both these parameters define, in some sense, the same thing, why not just club them together into a single parameter which gives the net spread?

The reason is 3 fold-

1. The **neighborhood is a function of the agents' behavior**, the **virulence is a function of the pathogen's behaviour**. It makes sense to keep them separate to deal with different cases of spread changes.
2. The model is such that there are **geometric differences** that arise by changing the neighborhood radius. The _buffering_ action of the refactory line is reduced by a larger neighborhood. It is not possible to increase the infection rate without also reducing this phenomenon if both were clubbed together.
3. The code would need too many changes.

!!! tips "Relevant quote"
    Any sufficiently advanced bug is indistinguishable from a feature.
    
    \- **Rich Kulawiec**

"""

# ╔═╡ 21b80dc0-67f9-11eb-0e87-6db59c88c296
md"""## Lockdowns and Neighborhoods

From the above discussion, it must be clear that a lockdown can be simulated pretty easily by varying the size of the neighborhood. The lockdown must also start at some predetermined point during the simulation. The size of the neighbourhood is changed by changing `grid.neighborhood` and the in-simulation change of the neighborhood can be executed by using the `postrunhook` of the `simulate!` function

"""

# ╔═╡ ead37ab2-66f9-11eb-2e56-07f105753991
function lockdown(grid, step; step_dict, radii)
    push!(radii, radius(grid.neighborhood))
    if step in keys(step_dict)
        grid.neighborhood = VonNeumannNeighborhood(step_dict[step], 2); 
    end
end

# ╔═╡ 5641a270-684b-11eb-29b4-ef8c59ef9679
md"## 200 → 1"

# ╔═╡ 0f20b710-6840-11eb-1712-9f4614e9ed87
ld1 = let
	step_dict = Dict(200 => 1) # Step -> radius
	radii = []
	grid = LifeGrid(
		(150, 150),
		4, 5,
		max=9,
		init=5,
		nh=VonNeumannNeighborhood(3, 2));
	results = simulate!(grid, 400; 
		postrunhook=lockdown, step_dict=step_dict, radii=radii);
	grid, results, radii
end;

# ╔═╡ 6f9dc770-6847-11eb-14ff-2b2e8562076f
@gif for frame ∈ ld1[2]
	plot(frame .|> p -> color(p, ld1[1]))
end every 2

# ╔═╡ 2bf1c112-6848-11eb-3b4a-a3ce569da696
makegraph(ld1...)

# ╔═╡ a9588a00-684b-11eb-2585-09a2cdf9ad80
md"## 50 → 2; 150 → 1; 300 → 3"

# ╔═╡ 4f4fdfd0-684c-11eb-3769-49d1677a1d52
ld2 = let
	step_dict = Dict(50 => 2, 150 => 1, 300 => 3) # Step -> radius
	radii = []
	grid = LifeGrid(
		(150, 150),
		4, 5,
		max=9,
		init=5,
		nh=VonNeumannNeighborhood(3, 2));
	results = simulate!(grid, 400; 
		postrunhook=lockdown, step_dict=step_dict, radii=radii);
	grid, results, radii
end;

# ╔═╡ d97dddfe-684c-11eb-1d11-5734830507a2
@gif for frame ∈ ld2[2]
	plot(frame .|> p -> color(p, ld2[1]))
end every 2

# ╔═╡ e1baa1c0-684c-11eb-292c-3336d6d715b1
makegraph(ld2...)

# ╔═╡ 52123eb0-684d-11eb-0179-895599ae7a2f
md"## 100 → 2; 300 → 3"

# ╔═╡ 4e182db0-684d-11eb-0bab-b7c22d3f3ff5
ld3 = let
	step_dict = Dict(100 => 2, 300 => 3) # Step -> radius
	radii = []
	grid = LifeGrid(
		(150, 150),
		4, 5,
		max=9,
		init=5,
		nh=VonNeumannNeighborhood(3, 2));
	results = simulate!(grid, 400; 
		postrunhook=lockdown, step_dict=step_dict, radii=radii);
	grid, results, radii
end;

# ╔═╡ 53bb3e50-684e-11eb-1bd3-d582200be8b5
@gif for frame ∈ ld3[2]
	plot(frame .|> p -> color(p, ld3[1]))
end every 2

# ╔═╡ 660d5660-684e-11eb-27ad-cbde3f616410
makegraph(ld3...)

# ╔═╡ f0c45ab2-684e-11eb-24ca-c5942b62f396
md"""## Observations

Apart from being very flashy, we see an interesting phenomenon. When the radius of the neighbourhood is changed, it takes about 25-30 steps for the grid to find a new equilibrium. During this transient phase, the distribution of states can vary wildly.

Further, the equilibrium fractions of the states averaged over a few cycles is nearly equal, but the exact variation may vary on the basis of initial conditions ([See 100 → 2; 300 → 3](#52123eb0-684d-11eb-0179-895599ae7a2f)). This can be attributed to the geometry of the neighborhood, the refactory buffer and the distribution of infected agents.

Because most large scale structures don't last for long times, we can assume that the final state distribution is memoryless or has only short-term memory. Thus our study of the effectiveness of a lockdown can be reduced to a simpler study of equilibrium infection as a function of neighborhood radius and [`grid.max`](#2563e210-67d5-11eb-089e-f376da973765), which, as you probably guessed it, is our next section.
"""

# ╔═╡ 61e28250-6857-11eb-39f9-2ff5565519ab
md"""
# Equilibrium infection distribution
**As a function of Neighbourhood Radius and Virulence of pathogen**
"""

# ╔═╡ 2734ec90-6859-11eb-2430-e3bd97302fd8
md"""## Introduction
This section is mainly computational. With the functionalities we have built up in the previous sections, we simply simulate a random 200X200 grid for 200 steps for different values of neighbourhood radius and virulence, find the average(over 5 cycle-lenghts) infected fraction then plot a scatter plot to see the distribution. 
"""

# ╔═╡ f7bf5830-6865-11eb-29d7-2be3d482ab20
function ss(grid, step; kwargs...)
	if count(i -> infectious(i, grid), state(grid)) == 0
		return :shortcurcuit
	end
end

# ╔═╡ 2e5c48b0-6859-11eb-2b69-0b78043bfc13
function equilibrium_infection(nhradius, p_max)
	arr = rand(1:9, 100, 100)
	grid = LifeGrid(arr, 4, 5, max=p_max, nh=VonNeumannNeighborhood(nhradius, 2))
	results = simulate!(grid, 200, postrunhook=ss)
	mean(count.(p -> infectious(p, grid), results[end - 45:end])) / 10000
end

# ╔═╡ 16c3dea0-685b-11eb-1b89-1da8ccf00941
md" ### Verification
We verify that the function is working correctly by checking for the case `nhradius = 3, p_max = 3` against the run [100 → 2; 300 → 3](#52123eb0-684d-11eb-0179-895599ae7a2f)"

# ╔═╡ ba6cd7a0-685b-11eb-3375-1b00d96dc710
begin
	res = equilibrium_infection(3, 9)
	if 0.35 < res < 0.4
		Markdown.parse("""
!!! tips "Success!"
    \$0.35< \\texttt{equilibrium\\_infection(3, 9)} = $(round(res, sigdigits=2)) <0.4\$
""")
	else
		Markdown.parse(md"""
!!! warning "Oh No!"
    \$0.35 \\nless \\texttt{equilibrium\\_infection(3, 9)} = $(round(res, sigdigits=2)) \\nless 0.4\$
    **This should never HAPPEN!!!**
""")
	end
end

# ╔═╡ 977a4062-6861-11eb-17eb-03874dfb3c79
md"""
!!! danger "Long running Code!"
    The following cell takes a long time (≈ 400 secs) to run. Be careful!
    You may want to reduce `nhs` and `maxs`, or skip some elements.
"""

# ╔═╡ 5b1e6e50-6877-11eb-1ca4-3df17b85575e
md"## Results"

# ╔═╡ 104a9570-685e-11eb-06c2-f9df4e63d1ba
begin
	nhs = 1:4
	maxs = 1:2:25
	points = [equilibrium_infection(nhradius, p_max) for nhradius in nhs, p_max in maxs]
end;

# ╔═╡ c6754200-685e-11eb-319f-f1d06650a620


# ╔═╡ 3f6db120-685e-11eb-3fa6-fd5ec5dfe981
let
	plotlyjs()
	p = plot(maxs, nhs, points, st=:surface, ylabel="NH Radius", xlabel="Virulence(max)", zlabel="Frac infected", title="Surface Plot")
	gr()
	p
end

# ╔═╡ 8af3d922-686e-11eb-03bb-d96c3ac24780
plot(maxs, nhs, points, ylabel="NH Radius", xlabel="Virulence(max)", title="Contour Plot", width=3)

# ╔═╡ 97bff6c2-6873-11eb-2f3d-c108ce4accde
begin
	plotlyjs()
	p = plot(maxs, nhs, points, ylabel="NH Radius", xlabel="Virulence(max)", title="HeatMap", st=:heatmap, hover = string.(reshape(points, :)))
	gr()
	p
end

# ╔═╡ 44881ce0-6877-11eb-0cef-f35b60baf357
md"
## Observations
1. As expected, runs with low virulence and low radius (bottom-right) have very low equilibrium infection levels.
2. For very virulent pathogens, even very strong lockdowns do not have much of an effect.
3. Even for mid-level virulence pathogens, weak lockdowns are only able to reduce the infection probability, but are unable to completely kill off the pathogen.
4. As we have seen before, even very low initial viral load can quickly explode in high interaction neighborhoods, so when the lockdown is removed, the infection quickly increases to pre-lockdown levels

## Discussion on the effectiveness of lockdowns

Given that lockdowns don't help in killing off the pathogen, was the whole last year a waste? Was it all just a government hoax to install 5G towers?

While the author can't say much about government plans to install 5G towers, it is clear that lockdowns have a dampening effect in the number of active cases. This decrease in cases allow the society to focus on combating the disease by spending more efforts into research for vaccines, treatments and preventive steps, along with building infrastructure and materials which are necessary to combat the recurrence of the pathogen. It also gives epidemologists and people working in Complex Systems time to make pretty plots and animations.
"

# ╔═╡ c7ea4fa0-687b-11eb-0dbc-7f316bd4aefd
md"# Further scope of this work
While this work has largely focused on the case of pathogen spread, the same work can be extended within the same topic and to other related topics.

## Things within the same topic
1. Fitting a mathematical model of the effect of lockdown and pathogen spread to data obtained by this work, then fitting it for other real life pathogens
2. Verification of model reflecting real pathogens via case studies
3. More extended statistical analysis to quantify equilibrium
4. Search for multi-equilibrium states
5. Contact tracing and quarantine simulations and comparing effectiveness of blanket lockdowns with more surgical methods.

## Other topics which can gain from this work
While we have called agents as people and the spread as a pathogen, any system which relaxes after some time can be modelled with a similar setup. Accordingly the meaning of lockdown can change.
Some examples include -
1. Spread of radioactivity with lockdown as shielding.
2. Complex web of processors which heat up after a while. There, the lockdown would be the tradeoff between a less complex web but a less efficient system of workload balancing.
3. Spread of algae and lockdown being the speed of flow of water.
"

# ╔═╡ 01866aa0-6854-11eb-3369-7deb6c18c3fe
md"""# References
[^Moitra2019]: Moitra, P., Sinha, S. Localized spatial distributions of disease phases yield long-term persistence of infection. Sci Rep 9, 20309 (2019). [https://doi.org/10.1038/s41598-019-56616-3](https://doi.org/10.1038/s41598-019-56616-3)
[^julia]: Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah. (2017) SIAM Review, 59: 65–98. doi: [10.1137/141000671](https://dx.doi.org/10.1137/141000671)
[^plots]: Plots - powerful convenience for visualization in Julia. Thomas Breloff. [docs.juliaplots.org](docs.juliaplots.org)
[^images]: JuliaImages: image processing and machine vision for Julia. [https://juliaimages.org/stable/](https://juliaimages.org/stable/)
[^pluto]: Pluto - Reactive Notebooks for Julia. Fons van der Plas and Mikołaj Bochenski. [https://github.com/fonsp/Pluto.jl](https://github.com/fonsp/Pluto.jl)
[^plutoui]: PlutoUI. Fons van der Plas. [https://github.com/fonsp/PlutoUI.jl](https://github.com/fonsp/PlutoUI.jl)
[^stats]: Statistics - Julia stdlib for Statistics. [https://docs.julialang.org/en/v1/stdlib/Statistics/](https://docs.julialang.org/en/v1/stdlib/Statistics/)
[^cellularbase]: CellularBase.jl: Abstract definitions for Cellular Automaton. Dhruva Sambrani. [https://dhruvasambrani.github.io/idc621/CellularBase](https://dhruvasambrani.github.io/idc621/CellularBase)
"""

# ╔═╡ Cell order:
# ╟─9a74ef90-67ad-11eb-1fc3-e7080fbfc7e1
# ╟─eb2aad32-687f-11eb-0c6b-91b443cdb3d3
# ╟─6f7bce82-67ad-11eb-1ec1-79ce2620fe9b
# ╟─3459c040-6858-11eb-09c1-95c7bdb02722
# ╟─7cb22830-67c9-11eb-0e0e-41f2316ed468
# ╟─65ff1d80-66f9-11eb-11c2-69c47b8458e4
# ╟─330799c0-6884-11eb-230d-71c512e4a447
# ╠═1bc961f0-66f7-11eb-17c8-2b54d4c4328d
# ╟─db4fd9e0-6883-11eb-384e-6df2c1d7a52b
# ╠═f8149ae0-66f6-11eb-0335-39baafdb24c8
# ╟─a95c8c70-67ad-11eb-124b-7b66d5547321
# ╟─73ea0ae0-66f9-11eb-19cf-b7ab6151ef63
# ╟─038cbdf2-6880-11eb-1b2b-2b4c4974be49
# ╠═35d0be40-66f7-11eb-2f60-01982265ccb0
# ╟─91afcec0-66f9-11eb-37b1-f78de5e0ba8a
# ╟─fcc23790-686e-11eb-03d7-1f0946446f49
# ╠═263db440-66f9-11eb-03b5-1fce784fad2d
# ╟─b8e7c920-66f9-11eb-397f-3760662a648c
# ╟─8b6833e0-6870-11eb-2593-69edaec0a7b8
# ╠═21d1eb20-66f8-11eb-26d3-4fb244ca8646
# ╟─de865ed2-66f9-11eb-31dc-c782d4338af3
# ╠═c3e6cc10-66f7-11eb-039f-df82addf1a5e
# ╠═4d4cdcc0-67b0-11eb-2278-8de64128d6f1
# ╠═7cdebc60-684b-11eb-38a2-45873e4b51d0
# ╠═5031733e-66f9-11eb-107f-6d69089d9384
# ╟─c0318a00-67ac-11eb-39b7-e3c13a5aff57
# ╟─47b81340-67cb-11eb-372e-972e4d4ff746
# ╠═02c24a70-67ae-11eb-2944-61a7963eb947
# ╟─fbdba862-67cf-11eb-39bc-f37eb3662785
# ╠═fe0dda70-66fa-11eb-0e02-73c05666aace
# ╠═8315d810-677a-11eb-3953-8faad42ada40
# ╟─7d24e300-67d0-11eb-0864-41c6166f13ac
# ╠═b8523010-67ce-11eb-2396-4780284dfa01
# ╠═cf22a450-67ce-11eb-2141-fd6f3efc689a
# ╟─07c440a0-67d1-11eb-3a8b-d559cac6276e
# ╠═30ce7a10-67bd-11eb-027d-8dc583201bc1
# ╠═4ca7fe50-67bd-11eb-382c-1304f3c1c3c4
# ╟─643ffb90-67ad-11eb-20c6-7744050e5448
# ╟─1cf29850-66fa-11eb-3f76-db755c0a4c50
# ╟─393a5bc0-67cb-11eb-3be1-a9ec7943b5f3
# ╟─2563e210-67d5-11eb-089e-f376da973765
# ╟─dafac300-67e9-11eb-0972-f9070193399e
# ╠═46055e90-67e9-11eb-2b6f-dffa58138005
# ╠═5626db00-67e9-11eb-0e12-a77b6b9c58be
# ╟─d4c2cff0-67e9-11eb-1f5b-03c7bab91782
# ╠═3056cdae-67ec-11eb-0b35-b3d6cde9b0fa
# ╠═3b0024a0-67ec-11eb-3398-eb07250758a9
# ╟─6d822d52-67ed-11eb-1b71-679d491bf133
# ╠═76e43730-67ed-11eb-201f-53bf8d3e3bfb
# ╠═8445d050-67ed-11eb-2411-f776a43c7f40
# ╟─0f5ba100-67ef-11eb-15f5-3fbaad70dde2
# ╟─9e71a7c0-67f1-11eb-0282-672d358d2292
# ╟─ef666b40-67f4-11eb-2f31-05145a6dfcba
# ╟─5586a260-67f4-11eb-11e2-b706fd0d4dfe
# ╟─21b80dc0-67f9-11eb-0e87-6db59c88c296
# ╠═ead37ab2-66f9-11eb-2e56-07f105753991
# ╟─5641a270-684b-11eb-29b4-ef8c59ef9679
# ╠═0f20b710-6840-11eb-1712-9f4614e9ed87
# ╠═6f9dc770-6847-11eb-14ff-2b2e8562076f
# ╠═2bf1c112-6848-11eb-3b4a-a3ce569da696
# ╟─a9588a00-684b-11eb-2585-09a2cdf9ad80
# ╠═4f4fdfd0-684c-11eb-3769-49d1677a1d52
# ╠═d97dddfe-684c-11eb-1d11-5734830507a2
# ╠═e1baa1c0-684c-11eb-292c-3336d6d715b1
# ╟─52123eb0-684d-11eb-0179-895599ae7a2f
# ╠═4e182db0-684d-11eb-0bab-b7c22d3f3ff5
# ╠═53bb3e50-684e-11eb-1bd3-d582200be8b5
# ╠═660d5660-684e-11eb-27ad-cbde3f616410
# ╟─f0c45ab2-684e-11eb-24ca-c5942b62f396
# ╟─61e28250-6857-11eb-39f9-2ff5565519ab
# ╟─2734ec90-6859-11eb-2430-e3bd97302fd8
# ╠═f7bf5830-6865-11eb-29d7-2be3d482ab20
# ╠═2e5c48b0-6859-11eb-2b69-0b78043bfc13
# ╟─16c3dea0-685b-11eb-1b89-1da8ccf00941
# ╟─ba6cd7a0-685b-11eb-3375-1b00d96dc710
# ╟─977a4062-6861-11eb-17eb-03874dfb3c79
# ╟─5b1e6e50-6877-11eb-1ca4-3df17b85575e
# ╠═104a9570-685e-11eb-06c2-f9df4e63d1ba
# ╟─c6754200-685e-11eb-319f-f1d06650a620
# ╠═3f6db120-685e-11eb-3fa6-fd5ec5dfe981
# ╠═8af3d922-686e-11eb-03bb-d96c3ac24780
# ╠═97bff6c2-6873-11eb-2f3d-c108ce4accde
# ╟─44881ce0-6877-11eb-0cef-f35b60baf357
# ╟─c7ea4fa0-687b-11eb-0dbc-7f316bd4aefd
# ╟─01866aa0-6854-11eb-3369-7deb6c18c3fe
