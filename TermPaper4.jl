### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 2a0680b2-a0f5-11eb-2658-d794c2f5d8a6
begin
	using LightGraphs
	using LinearAlgebra
	using Plots
	using PlutoUI
	using GraphRecipes
	using Statistics
	using LsqFit
	using Random
end

# ╔═╡ 2ea044d4-0ef4-4777-be1d-120b6e33b853
html"""<center><h1 style="font-size: 3.5rem">Term Paper 4<h1></center><br>
This work is a submission for the course Modelling Complex Systems(IDC-631) for the year 2020-21 Spring (Even) Semester by <b>Dhruva Sambrani, MS18163</b>
"""

# ╔═╡ 521ee215-35b3-432f-b584-1cc98daa4c5d
md"""
# Abstract

This paper looks at the robustness of the mean path length of many network topologies to removal of edges. A definition of "fitness" of a graph is also provided which attempts to capture real-world constraints on the networks we are considering.

## Structure of the paper

- **§1** Graphs and Networks - A brief introduction
- **§2** The Storm Story - The problem statement
- **§3** Gyms and Strength - A definition of network robustness
- **§4** Results and Conclusions
- **§5** References

!!! note "Note to reader"
    The author believes that papers which involve a lot of simulation and code should also be accompanied with the code they use so that readers can gain deeper insight into the processes described in the paper. Hence, the author has attempted to create a paper which interlaces the information with the corresponding code. However, in the interest of keeping the code aligned to the topic at hand, certain auxillary functions have been hidden. The content of these hidden cells can be found in the executable version of this notebook [here]()

The coding language used is `Julia`[^julia] and following document is made with `Pluto.jl` [^pluto].
"""

# ╔═╡ 5d01718d-aa67-467d-9a9c-1f65ac6d6dca
TableOfContents(depth=6)

# ╔═╡ 6b4e7cd4-d8ab-4d1b-a0df-a299595e4739
md"
# §1 Graphs and Networks - A brief introduction

Most things in our world do not exist in isolation, and are related to each other in some way or another. Take all people. They can be connected by, well, relationships. Houses are connected over a powergrid, and computers over a network. These networks are modelled mathematically by a structure known as Graphs. 

Graphs are completely definied by a set of vertices(nodes) $\{v_i\}$ and edges $\{e_{i\rightarrow j}\}$. Often, the vertices are considered to be individuals in the system and the edges are considered to be the existance of a relationship between the edges.
"

# ╔═╡ 7fd589b6-e188-4f67-a6eb-e04d4c886c51
md"""## Example Graphs
$(plot(
	graphplot(barbell_graph(5, 5), title="Barbell"), 
	graphplot(complete_graph(10), title="Complete"),
	graphplot(watts_strogatz(10, 4, 0.3), title="Watts-Strogtaz"),
	graphplot(erdos_renyi(10, 0.4), title="Erdos-Renyi"),
))"""

# ╔═╡ 1d8f1d47-d0f7-4901-b91f-4308f79681f5
md"""## Graph Properties

Graphs, like all other mathematical objects have certain properties associated with them.
"""

# ╔═╡ 63aed5e9-f17a-4468-839b-e5f01d39e334
md"""
### Adjacency Matrix

An adjacency matrix is a bit matrix which is defined as $A_{ij} = e_{i\rightarrow j}$. For undirected, unweighted graphs, this is a symmetric matrix. The adjacency matrix has all the information about the topology of the graph.

For example, the adjacency matrix of a Barbell graph looks like this-
"""

# ╔═╡ 5e5460e6-180f-4474-895a-8eee1e7c181a
adjacency_matrix(barbell_graph(5,5))

# ╔═╡ 92b62f06-c57d-43f5-87d7-a5ae409e1402
md"### Connectedness

A graph is connected if there exists a path between any two nodes on the graph.

As an example, an Erdos-Renyi graph may not be connected
"

# ╔═╡ 8aec44e0-3aa9-4228-87a4-d605d23e2ded
g = erdos_renyi(10,20)

# ╔═╡ 48a7f6cd-a457-4c3a-b47a-2ccc612b0878
graphplot(g, size=(300,225))

# ╔═╡ 069d1c41-f378-4364-ad1e-cd80c3dcd8e4
is_connected(g) ? md"**Graph is connected**" : md"**Graph is not connected**"

# ╔═╡ ff945901-8c13-4404-b71c-33616a478aef
md"""# §2 The Storm Story - The problem statement

Consider this. It is a stormy night when your workplace network stops working and you are cutoff from the world. You are only scrolling Social Media, and you imagine remember the good times when you could "like" your friends' pictures. Living in the 21st century, you'd expect the network to have multiple paths to the internet.

This paper asks exactly that. How many connections can be dropped from the network before some devices are completely disconnected from the rest of the world. While you could wait for another 10 minutes before your site loading, there are some situations which demand that the network stay up for as long as possible before complete failure. These include banking and security systems, emergency response networks, remote control systems and food ordering systems. Ok, the last one was a stretch.

Similar models can also be used to model divergence of species in a population where some geo-physical/biological barrier prevents genes from spreading between the sub-populations.

Along the same lines, power grids too can be tested for robustness. However, powergrids add an extra layer of complexity, as 
1. They are directed graphs
2. Source and sink nodes are not equivalent
Due to these complexities, we do not look at this particular problem and leave it for future work.

However, we can add certain further complexity to the problem. While the network may remain connected, the mean path length can increase to such high values, that it becomes _effectively_ disconnected. This is especially important in the case where there are competing processes at play. 
- In the case of genes, it is the rate of mutations and the rate of change of the environment.
- In the case of emergency response, it is the priority of the case
- In remote control, too high ping rates can lead to poor control
- In the case of internet, the competing process is often the patience of the user.

Hence, we need a way to measure the **robustness** of the network to connection failures.
"""

# ╔═╡ 89fdc6e6-a16a-4c53-b64f-409798d1be7f
md"""
# §3 Gyms and Strength - A definition of network robustness

## The Graph theory model

Graph theory is an elegant way to translate the problem into a mathematical and computational struture that can be manipulated.

In this scenario, the nodes are the devices, organisms, or some equivalent entity. The edges are the "connectivity" between two nodes. In the case of an internet network, it is the channel for information to pass through."""

# ╔═╡ 926e2dcd-5ee7-4507-b41d-879b58f42662
md"""## Network failure

As described in the previous section, the network "fails" when it becomes disconnected or mean path length increases beyond a particular cutoff.

!!! note "Network Failure"
	- A network **fails** when,
	  - It becomes disconnected OR
	  - Its mean path length increases beyond some cutoff, as a multiple of initial mean path length or an absolute value."""

# ╔═╡ ffc21b99-e466-424c-918a-d2cd30844278
md"## Topology dependance

Naively, it seems like that the network which has more edges will fare better than one which has few. But this is not the case. Consider the Barbell graph.
"

# ╔═╡ 5a5e98c6-d983-4fbb-b073-c347968a6415
graphplot(barbell_graph(5,5), size=(300,225))

# ╔═╡ e4c5c26b-4976-44af-8aa9-2d7d22f84af0
md"In this graph, if the edge which connects the two subgraphs(cliques) is broken, the graph immediatly fails"

# ╔═╡ 01455f58-0492-41a6-ae44-3cd39750aaa1
let
	g = barbell_graph(5,5)
	rem_edge!(g, 5, 6)
	isc = is_connected(g) ? "connected" : "not connected"
	md"$(graphplot(g, size=(300,225))) Graph is $isc"
end

# ╔═╡ b19b1856-c46b-4812-bc3e-9e7dce4b5b69
md"Hence, the robustness of the graph is most definitely dependant on the topology of the network."

# ╔═╡ 8fe47f6d-1a41-4e58-b7dd-eb365065404f
md"""## Robustness

The simulation then could involve creating a graph, checking connectedness, breaking edges one by one and checking connectedness after each edge removal. Then, a **measure of robustness is the amount of edges we can remove before the network fails**. However, it would be unfair to a graph of 100 edges which fails on removing 50 edges to be considered worse than a network of 1000 edges that fails on removing 100 edges. So a better measure of robustness is the **_fraction_ of edges which can fail before the network fails**. As we saw in the previous subsection, robustness clearly depends _which_ edge we remove. Hence, we again amend our definition of robustness to read as **_the average_ fraction of edges which can fail before the network fails**

!!! note "Definition of Robustness"
	Robustness of a network is defined as the the _average fraction_ of edges which can fail before the network _fails_.
"""

# ╔═╡ 07953ab8-6449-4968-8d3d-75548706170e
md"""### Small scale Robustness

Sometimes, it may be more useful to measure how the mean path length increases if few edges are removed, especially when edge failure is a rare event, or faults can be fixed quickly, and when responsiveness to the network is of higher importance than long scale robustness. For this, we define small scale robustness of a network.

!!! note "Small scale robustness"
	It is defined as the slope of the mean path length vs the fraction of nodes removed for small fractions.

There are some evident disadvantages to such a definition
- Definition of "small" - this can be solved by choosing a range of interest
- Difficult to define for small networks - Since the definition involves fitting, small networks suffer having too few data points which can not be increased.
- Computationally expensive - Compared to our definition of robustness, SSR is more computationally expensive.
- Linear nature of plot is not guaranteed - Consider the Barbell network. 
- Linear nature may not remain over the range of interest - More generalized fits may be more useful
"""

# ╔═╡ 0a470cac-b6a0-45bf-9b67-b1ad22458011
md"# §4 Results and Conclusion

In the following section, we implement the setup above for multiple topologies, and look at which topologies are more robust according to our definition.
"

# ╔═╡ 611f678e-48a7-4a0b-8f8d-9abfea658eba
md"## Implementation"

# ╔═╡ 894715a6-9ca5-49fb-be98-ad9131efbbfc
begin
	# Remove a random edge
	breakrandom!(g::AbstractGraph) = rem_edge!(g, rand(collect(edges(g))))
	mean_path_length(g) = 
		sum(floyd_warshall_shortest_paths(g).dists)/(nv(g)^2-nv(g))
end

# ╔═╡ 2902dc16-2a73-42bd-b2bb-7e43e3803c39
md"### Mean Path Length

This is the mean of the shortest path lengths between all nodes in the graph.

It is given by the formula 

$$\mu(g) = \frac{\sum_{i\in\text{v}(g), j\in\text{v}(g),\ j\ne i} d(i,j)}{n(n-1)}$$

where $d(i,j)$ is the shortest possible path $i\rightarrow j$. If no such path exists, some definitions make it infinity, others make it 0. The appropriate definition should be chosen as per requirements.

The Floyd Warshall shortest paths algorithm provides a matrix of the shortest paths between all the nodes in a graph.

For example, the mean path length of a completely connected graph is $(mean_path_length(complete_graph(10))), as expected.
"

# ╔═╡ bef640f9-95c9-4e74-b7ad-91ac793f3a97
function ssr(pathlengths, frac)
	datarange = floor(Int, length(pathlengths)*frac)
	datarange > 5 || return -1 # network robustness is too low to measure SSR
	@. model(x, (m, c)) = m*x+c
	fit = curve_fit(
		model, 
		collect(range(0, stop=frac, length=datarange)), 
		pathlengths[1:datarange], 
		[0.,0.]
	)
	return fit.param[1]
end

# ╔═╡ 8e9b6375-bfc4-4e55-9afd-a6e8c810b900
begin
	mpl_increase(mpls) = last(mpls)/first(mpls)
	function simulate!(g, cutoff_mpl_scale; seed=-1)
		if seed >= 0
			Random.seed!(seed)
		end
		pathlengths = [mean_path_length(g)]
		nedges = ne(g)
		while is_connected(g) || (mpl_increase(pathlengths) > cutoff_mpl_scale)
			breakrandom!(g)
			push!(pathlengths, mean_path_length(g))
		end
		pathlengths = length(pathlengths)>1 ? pathlengths[1:end-1] : pathlengths
		return (
			length(pathlengths)/nedges, 
			pathlengths, 
			length(pathlengths)>1 ? range(0, stop=length(pathlengths)/nedges, length=length(pathlengths)) : [0]
		)
	end
end

# ╔═╡ c71ed1c4-a13b-432c-93a3-c9d27cc6ce94
function makeplot((pls, x), name)
	plot(
		x,
		map(pls) do i;
			i/first(pls);
		end,
		xlims=(0,1),
		ylims=(1,5),
		xlabel="Fraction of edges removed",
		ylabel="Ratio of MPL",
		title=name
	)
end

# ╔═╡ c4e14fd7-d0de-4370-ad7f-06563cff1009
md"## Test Runs"

# ╔═╡ 22c6f858-d074-49ae-bf41-b3502c08ed15
seed = 5176

# ╔═╡ a1d6a627-fc7c-4e9c-9a64-14c6059dbb3f
begin
	p1 = makeplot(simulate!(watts_strogatz(50, 6, 0.1, seed=seed), Inf, seed=seed)[2:3], "Watts-Strogatz")
	p2 = makeplot(simulate!(complete_graph(50), Inf, seed=seed)[2:3], "Complete")
	p3 = makeplot(simulate!(barabasi_albert(50, 6, seed=seed), Inf, seed=seed)[2:3], "Barabasi-Albert")
	plot(p1, p2, p3, size=(600, 600))
end

# ╔═╡ fda41875-dced-4ea3-953a-b0691fc3aa63
function get_robustness_measures(generator, steps; cutoff_mpl_scale=Inf, frac=0.1)
	rbms = map(1:steps) do i
		g = generator()
		data = simulate!(g, cutoff_mpl_scale)
		println(length(data[2]))
		(data[1], ssr(data[2], frac))
	end
	return first.(rbms), filter(i->i>0, last.(rbms))
end

# ╔═╡ 57b518a2-be25-47e6-b556-f1707c35405b
begin
	function getstats(rbs, ssrs, name)
		ssrs = length(ssrs)==0 ? [-1] : ssrs
		md"""
	### $(name)
	
	| Property | Value for Robustness             |Value for SSR                      |
	|:--------:|:--------------------------------:|:---------------------------------:|
	| Median   |$(round(median(rbs), sigdigits=3))|$(round(median(ssrs), sigdigits=3))|
	| Mean     |$(round(mean(rbs), sigdigits=3))  |$(round(mean(ssrs), sigdigits=3))  |
	| StdDev   |$(round(std(rbs), sigdigits=3))   |$(round(std(ssrs), sigdigits=3))   |
	| Length   |$(length(rbs))                    |$(length(ssrs))                    |
	"""
	end
	getstats(name) = (dat) -> getstats(dat[1], dat[2], name)
end

# ╔═╡ 6370353a-fbc4-4f1b-be96-1942b486e617
md"## Runs over different topologies"

# ╔═╡ df107cca-e86f-48a8-9c7a-6c9c6036bc94
get_robustness_measures(()->watts_strogatz(50, 6, 0.01), 500) |> getstats("Watts-Strogatz low β(0.01)")

# ╔═╡ b1468ecd-c3e8-4009-9290-e944bb433ef5
get_robustness_measures(()->watts_strogatz(50, 6, 0.1), 500) |> getstats("Watts-Strogatz mid β(0.1)")

# ╔═╡ 8b0af9d6-e60b-483e-a848-0a1d8963e351
get_robustness_measures(()->watts_strogatz(50, 6, 0.5), 500) |> getstats("Watts-Strogatz high β(0.5)")

# ╔═╡ 576d7b10-8fbf-4aa8-893e-5dbf7b567812
get_robustness_measures(()->watts_strogatz(50, 20, 0.1), 500) |> getstats("Watts-Strogatz high k(20)")

# ╔═╡ 3ac2fc49-b17e-48d1-8991-d7f970331c31
get_robustness_measures(()->barabasi_albert(50, 6), 500) |> getstats("Barabasi-Albert mid k(6)")

# ╔═╡ 64b84265-b49d-4869-adcb-b9c68d505c4f
get_robustness_measures(()->barabasi_albert(50, 20), 500) |> getstats("Barabasi-Albert high k(20)")

# ╔═╡ 58be2740-f3ab-46af-9d64-b3897a7deeb3
get_robustness_measures(()->complete_graph(50), 100) |> getstats("Complete Graph")

# ╔═╡ 51a0b307-cb71-4280-af8f-8ad5afafa8c9
get_robustness_measures(()->barbell_graph(25,25), 100) |> getstats("Barbell Graph")

# ╔═╡ 4c13a4fe-f642-469b-8203-0d7796457581
get_robustness_measures(()->erdos_renyi(50, 100), 100) |> getstats("Erdos-Renyi Graph low k(10%)")

# ╔═╡ f35f13d3-03d3-4ed1-9246-4ae0f06835ac
get_robustness_measures(()->erdos_renyi(50, 625), 100) |> getstats("Erdos-Renyi Graph high k(50%)")

# ╔═╡ 67c579be-f039-4352-bb9b-1791ec1ec858
md"""## Conclusions

1. As a comparision between kinds of topology, a Complete graph has the highest robustness. This is not surprising, because after all, there are so many connections that random edges being removed don't affect the topology of the graph, and no particular edge increases leads to network failure.
2. For the Erdos-Renyi graph[^erdos_renyi], a very low k can lead to multiple nodes outside the main graph. For a high k results are similar to that of the Complete graph, which is not surprising, as a complete graph with random nodes removed is a random graph.
3. For the Watts-Strogatz[^watts_strogatz] graph, while increasing β decreases SSR, it also decreases robustness. If k is increased, it increases robustness and decreases SSR. This is in contrast to simply adding more nodes in the Erdos-Renyi, as it decreases SSR, which is preferable. This strongly shows that Watts-Strogatz model is resistant to small scale perturbations, leading to an evolutionary benefit, as predicted before by [^wiki].
4. The Barabasi-Albert[^barabasi_albert] graph is interesting in context of the other graphs, especially when compared to the Watts-Strogatz model. While Barabasi-Albert is more robust and has lower SSR for low k, for higher ks, the Barabasi-Albert has lower robustness and comparable SSR. This suggests that for low density edge graphs, a Barabasi-Albert topology is "better" than the Watts-Strogatz, but as density increases, the Watts-Strogatz model scales better.
5. The Barbell Graph reacts surprisingly well, but it has low robustness and high SSR.
"""

# ╔═╡ f7393bf5-ad86-4619-a69a-9ff7108b2870
md"# References
[^robustness]:[wiki: Small World Network#Network_robustness](https://en.wikipedia.org/wiki/Small-world_network#Network_robustness)
[^erdos_renyi]: [wiki: Erdos-Renyi_model](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model)
[^watts_strogatz]:[Watts-Strogatz_model](https://en.wikipedia.org/wiki/Watts%E2%80%93Strogatz_model)
[^barabasi_albert]: [wiki: Barabasi-Albert_model](https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model)

## Packages used 

[^lightgraphs]: JuliaGraphs/LightGraphs.jl: an optimized graphs package for the Julia programming language. Seth Bromberger, James Fairbanks, and other contributors. [10.5281/zenodo.889971](https://doi.org/10.5281/zenodo.889971)
[^graphrecepies]: GraphRecipes. Thomas Breloff, and other contributors. [github](https://github.com/JuliaPlots/GraphRecipes.jl)
[^julia]: Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah. (2017) SIAM Review, 59: 65–98. doi: [10.1137/141000671](https://dx.doi.org/10.1137/141000671)
[^plots]: Plots - powerful convenience for visualization in Julia. Thomas Breloff. [docs.juliaplots.org](docs.juliaplots.org)
[^pluto]: Pluto - Reactive Notebooks for Julia. Fons van der Plas and Mikołaj Bochenski. [https://github.com/fonsp/Pluto.jl](https://github.com/fonsp/Pluto.jl)
[^plutoui]: PlutoUI. Fons van der Plas. [https://github.com/fonsp/PlutoUI.jl](https://github.com/fonsp/PlutoUI.jl)
[^stats]: Statistics - Julia stdlib for Statistics. [https://docs.julialang.org/en/v1/stdlib/Statistics/](https://docs.julialang.org/en/v1/stdlib/Statistics/)
[^lsqfit]: LsqFit.jl: Basic least-squares fitting in pure Julia under an MIT license. [https://github.com/JuliaNLSolvers/LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl)
"

# ╔═╡ Cell order:
# ╟─2ea044d4-0ef4-4777-be1d-120b6e33b853
# ╟─521ee215-35b3-432f-b584-1cc98daa4c5d
# ╟─5d01718d-aa67-467d-9a9c-1f65ac6d6dca
# ╟─2a0680b2-a0f5-11eb-2658-d794c2f5d8a6
# ╟─6b4e7cd4-d8ab-4d1b-a0df-a299595e4739
# ╟─7fd589b6-e188-4f67-a6eb-e04d4c886c51
# ╟─1d8f1d47-d0f7-4901-b91f-4308f79681f5
# ╟─63aed5e9-f17a-4468-839b-e5f01d39e334
# ╟─5e5460e6-180f-4474-895a-8eee1e7c181a
# ╟─2902dc16-2a73-42bd-b2bb-7e43e3803c39
# ╟─92b62f06-c57d-43f5-87d7-a5ae409e1402
# ╟─8aec44e0-3aa9-4228-87a4-d605d23e2ded
# ╟─48a7f6cd-a457-4c3a-b47a-2ccc612b0878
# ╟─069d1c41-f378-4364-ad1e-cd80c3dcd8e4
# ╟─ff945901-8c13-4404-b71c-33616a478aef
# ╟─89fdc6e6-a16a-4c53-b64f-409798d1be7f
# ╟─926e2dcd-5ee7-4507-b41d-879b58f42662
# ╟─ffc21b99-e466-424c-918a-d2cd30844278
# ╟─5a5e98c6-d983-4fbb-b073-c347968a6415
# ╟─e4c5c26b-4976-44af-8aa9-2d7d22f84af0
# ╟─01455f58-0492-41a6-ae44-3cd39750aaa1
# ╟─b19b1856-c46b-4812-bc3e-9e7dce4b5b69
# ╟─8fe47f6d-1a41-4e58-b7dd-eb365065404f
# ╟─07953ab8-6449-4968-8d3d-75548706170e
# ╟─0a470cac-b6a0-45bf-9b67-b1ad22458011
# ╟─611f678e-48a7-4a0b-8f8d-9abfea658eba
# ╠═894715a6-9ca5-49fb-be98-ad9131efbbfc
# ╠═bef640f9-95c9-4e74-b7ad-91ac793f3a97
# ╠═8e9b6375-bfc4-4e55-9afd-a6e8c810b900
# ╠═c71ed1c4-a13b-432c-93a3-c9d27cc6ce94
# ╟─c4e14fd7-d0de-4370-ad7f-06563cff1009
# ╠═22c6f858-d074-49ae-bf41-b3502c08ed15
# ╠═a1d6a627-fc7c-4e9c-9a64-14c6059dbb3f
# ╠═fda41875-dced-4ea3-953a-b0691fc3aa63
# ╟─57b518a2-be25-47e6-b556-f1707c35405b
# ╟─6370353a-fbc4-4f1b-be96-1942b486e617
# ╟─df107cca-e86f-48a8-9c7a-6c9c6036bc94
# ╟─b1468ecd-c3e8-4009-9290-e944bb433ef5
# ╟─8b0af9d6-e60b-483e-a848-0a1d8963e351
# ╟─576d7b10-8fbf-4aa8-893e-5dbf7b567812
# ╟─3ac2fc49-b17e-48d1-8991-d7f970331c31
# ╟─64b84265-b49d-4869-adcb-b9c68d505c4f
# ╟─58be2740-f3ab-46af-9d64-b3897a7deeb3
# ╟─51a0b307-cb71-4280-af8f-8ad5afafa8c9
# ╟─4c13a4fe-f642-469b-8203-0d7796457581
# ╟─f35f13d3-03d3-4ed1-9246-4ae0f06835ac
# ╟─67c579be-f039-4352-bb9b-1791ec1ec858
# ╟─f7393bf5-ad86-4619-a69a-9ff7108b2870
