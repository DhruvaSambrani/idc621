### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 8d7c0eb0-84ca-11eb-0018-bdb2e389ee52
begin
	using PlutoUI
	using Distributions
	using Statistics
	using Plots
	using Random
	using JLD
	using Images
	using LsqFit
	md"_hidden code cell_"
end

# ╔═╡ e8d31f7c-88aa-11eb-05b7-31ff25d50afb
html"""<center><h1 style="font-size: 3.5rem">Term Paper 3<h1></center><br>
This work is a submission for the course Modelling Complex Systems(IDC-631) for the year 2020-21 Spring (Even) Semester by <b>Dhruva Sambrani, MS18163</b>
"""

# ╔═╡ 0d6813c6-88ae-11eb-3f39-b59ffd71ffd5
md"# Abstract

This paper mainly attempts to recreate Kuramoto's model of coupled oscillators as defined in [^strogatz].
This paper consists of four sections- 

- §1. Describe and implement the model as in the paper
- §2. Check that the model evolves as we expect it to and make pretty animations
- §3. Show that simulated results agree with certain analytical properties
- §4. Study the robustness of the model to different distributions

"

# ╔═╡ 831e6f04-8a55-11eb-35da-2f77e5dd68a5
md"## How to read the paper

This paper is made using the Pluto[^pluto] Package which provides an interface which can show both code and theory side by side, which will hopefully enhance the understanding of the work and also provide a reproducable implementation without having to go into the details of it. Those who wish to read the paper for the theory and results can easily skip over the code and still understand everything. Those who are attempting to recreate the work can use the code as a guide. However, in the interest of keeping the code aligned to the topic at hand, certain auxillary functions have been hidden. The content of these hidden cells can be found in the executable version of this notebook [here](https://github.com/DhruvaSambrani/idc621/blob/main/TermPaper3.jl)
"

# ╔═╡ bd9326c2-8954-11eb-003e-b52f7107f199
TableOfContents()

# ╔═╡ 3effc012-8a56-11eb-1cd2-475e917c969b
md"# §0 Packages used

The Plots.jl[^plots] package is used for plotting and the Statistics[^stats] Package for some statistical analysis. Distributions.jl[^distributions] is used to sample the . PlutoUI.jl[^plutoui] is used for some UI elements such as the ToC. JLD[^jld] is used for storing data between runs. Images.jl[^images] is used for Gaussian smoothing. LsqFit.jl[^lsqfit] is used to fit the $\kappa_c$ vs $\sigma$ function."


# ╔═╡ 27e115ae-88b0-11eb-1d7d-57ae925c4c9f
md"# §1 - The Kuramoto Model

This is a model which attempts to simulate coupled oscillators. Consider a set of $n$ oscillators $A_{i\ |\ 1\le i\le n}$. Each $A_i$ has an angular position and internal angular velocity $\theta_i$ and $\omega_i$. Then, the most general coupled oscillator equation would be 

$$\dot{\theta}_j = \omega_j + f(\theta_i\dots)_{1\le i\le n}$$

If the coupling is weak, and symmetric, this can be reduced to 

$$\dot{\theta}_j = \omega_j + \kappa \sum_{i=1}^{n} \sin(\theta_j-\theta_i)$$

Which can be recast into the following equation

$$\dot{\theta}_j = \omega_j + \kappa r \sin(\Psi-\theta_i)$$
$$\text{where}$$
$$re^{\iota\Psi} = \text{mean}(e^{\iota\theta_i}\dots)_{1\le i \le n}$$
"

# ╔═╡ bc07ce24-88b5-11eb-2052-53543e8ef9fc
md"## Agent

Let us define a `struct`ure called Agent which holds two variables R and V, which correspond to $\theta$ and $\omega$ respectively.
"

# ╔═╡ a34a7c40-84ca-11eb-3433-41615dfc1c08
mutable struct Agent
	R::Float64
	V::Float64
end

# ╔═╡ 56ce942c-88b8-11eb-08fc-0970f8de6359
md"Now let's define how an agent updates, given κ, r, Ψ and dt, as defined in the equations before"

# ╔═╡ c4b72810-84ca-11eb-1855-23703a18b3e4
function update!(agent::Agent, κr::Float64, Phi::Float64, dt::Float64)
	agent.R += (agent.V + κr * sin(Phi - agent.R)) * dt
	# θ = θ + dθ/dt * Δt = θ + (ω + κr⋅sin(Ψ-θ))⋅dt
end

# ╔═╡ 9a1dec4e-88e7-11eb-0736-6b6e6b24627f
md"## Helper Functions"

# ╔═╡ c12bde70-84ca-11eb-1fe2-634e16d69909
begin
	Statistics.mean(arr::Array{Agent})::Complex = mean(map(arr) do ag; cis(ag.R); end)
	init(n::Int, σ::Float64) = Agent.(1., rand(Normal(0., σ), n))
end

# ╔═╡ 9fa5bf26-89b8-11eb-0636-ed01bb9b689b
init(n::Int, d::Distribution) = Agent.(1., rand(d, n))

# ╔═╡ a7d91872-88e7-11eb-2799-3502fc01e720
md"## Evolution and Simulation

As always, we define a single step evolution along with a multistep simulate function.

The evolve function takes a list of agents and updates them according to the dynamics explained above.
The simulate function sets the parameters and records certain variables for plotting and further analysis
" 

# ╔═╡ cb3e73a0-84ca-11eb-3c68-8fc3fc5cf239
function evolve!(agents::Array{Agent}, κ::Float64, dt::Float64)
	meanvec = mean(agents) 
	update!.(agents, κ*abs(meanvec), angle(meanvec), dt)
	meanvec
end

# ╔═╡ d00b64fe-84ca-11eb-1e68-f73bbca9e373
begin
	function simulate(steps::Int, n::Int, σ::Float64, κ::Float64, time::Int64)
		dt = time/steps
		agents = init(n, σ)
		means = Array{Complex,1}(undef, Int(steps))
		for step in 1:steps
			means[step] = evolve!(agents, κ, dt)
		end
		means
	end
	function simulate(steps::Int, agents::Array{Agent,1}, κ::Float64, time::Int64)
		dt = time/steps
		agents = deepcopy(agents)
		means = Array{Complex,1}(undef, Int(steps))
		for step in 1:steps
			means[step] = evolve!(agents, κ, dt)
		end
		means
	end
end

# ╔═╡ 393ed800-88e9-11eb-2160-edef4c90ee74
md"# §2 Visualizing the model

Let us first make sure that the model is running correctly. We can either visualize the positions of the oscillators with time, or we can plot $r$ vs time and make sure that the results match what is given in the paper. We will do both.
"

# ╔═╡ 86408954-88f9-11eb-2ee7-85532daf659d
md"## $r_\infty$ depends on κ"

# ╔═╡ 8c9b5c8e-84d3-11eb-0f9d-47c844a46a68
begin
	no_corr_means = simulate(100_000, 1000, 3., 1., 40)
	corr_means = simulate(100_000, 1000, 3., 7., 40)
end;

# ╔═╡ c618ba40-84d2-11eb-3990-93a786c8f8a9
begin
	plot(abs.(corr_means), ylims=(0,1), label="κ=7")
	plot!(abs.(no_corr_means), ylims=(0,1), label="κ=1")
	plot!(legend=:right, xlabel="Timestep", ylabel="r", title="r vs t")
end

# ╔═╡ cbb28294-88f9-11eb-1b77-fdacf61b8c9b
md"## Animated visualization of the dynamics

For this, we define a `simulate2` function which is of the same form as the previous function, but it also records the agents for each time step. 
Then we plot the locations of the oscillators and the mean vector for each time step and make an animation.
Due to computational limitations, only 100 particles for 10000 timesteps are taken, and animation is made for every 20th frame. Also, for reproducability, the random seed is fixed. The seed was selected so that the results were clearer to show. However all observations are observed for any number of unseeded runs. To test, simply comment out the seed line.
"

# ╔═╡ f564ad60-88ea-11eb-31c5-a9f36ba1c93b
function simulate2(steps::Int, n::Int, σ::Float64, κ::Float64, time::Int64)
	dt = time/steps
	agents = init(n, σ)
	means = Array{Complex,1}(undef, Int(steps))
	recd = Array{Array{Agent, 1}, 1}(undef, Int(steps))
	for step in 1:steps
		means[step] = evolve!(agents, κ, dt)
		recd[step] = deepcopy(agents)
	end
	means, recd
end

# ╔═╡ bef4689a-88ec-11eb-2297-3bfe711d504b
begin
	Random.seed!(1901)
	ms, ags = simulate2(10000, 100, 3., 7., 40)
end;

# ╔═╡ 2a5fcef8-88ed-11eb-36e9-a743d38881f6
begin
	cmpl = ags .|> a->map(a) do ag
		cis(ag.R)
	end
	anim = @animate for i in 1:20:10000
		plot(real.(cmpl[i]), imag.(cmpl[i]), st=:scatter, xlims=(-1,1), ylims=(-1,1), aspect_ratio=1)
		plot!([0,real(ms[i])], [0, imag(ms[i])], width=2, label="Mean Vector")
	end
	md"_hidden code cell_"
end

# ╔═╡ 586e0782-88f3-11eb-2f30-03211c7a70a2
gif(anim, fps=3)

# ╔═╡ f2a5c968-88fb-11eb-3a39-2b95d2a8d978
begin
	plot(abs.(ms), ylim=(0,1), label="r", title="r vs Timestep", xlabel="Timestep", ylabel="r")
end

# ╔═╡ 222d4b00-88f4-11eb-2f5b-bd320a59b2ab
begin
	rs = ags .|> a->map(a) do ag; ag.R; end
	vs = (rs[2:end] .- rs[1:end-1]) .* 250
	md"_hidden code cell_"
end

# ╔═╡ 9537f2f8-88f4-11eb-091d-b5827649cb40
plot(1:10:10000, hcat(vs...)'[1:10:10000, :], legend=false, xaxis="Timestep", yaxis="Velocity", title="Effective Velocity vs time")

# ╔═╡ def3ccb4-8949-11eb-1ab1-a3403dc946ad
md"# §3 Simulations agree with analysis

In section 3.3 of [^strogatz], the author mentions the following observations of the simulations."

# ╔═╡ d2d048a8-8956-11eb-0536-972068f6862f
md"## r∞ vs κ

- For all $K$ less than a certain threshold $K_c$, the oscillators act as if they were uncoupled: the phases become uniformly distributed around the circle, starting from any initial condition. Then $r(t)$ decays to a tiny jitter of size $\mathcal O(N^{−1/2})$
- When $K > K_c$, this incoherent state becomes unstable and $r(t)$ grows exponentially, reflecting the nucleation of a small cluster of oscillators that are mutually synchronized, thereby generating a collective oscillation. Eventually $r(t)$ saturates at some level $r_\infty < 1$, though still with $\mathcal O(N^{−1/2})$ fluctuations.

Both of these claims are shown to be true in the first plot. For a \"high\" $K(=7)$, $r$ settles with jitter at 1, and for a \"low\" $K(=1)$, $r$ settles with jitter at 0.
"

# ╔═╡ 963c13ba-8957-11eb-3d45-f1d14f7b1c87
md"
## Individual Oscillators

- The population splits into two groups: the oscillators near the center of the frequency distribution lock together at the mean frequency and co-rotate with the average phase ψ(t), those in the tails run near their natural frequencies and drift relative to the synchronized cluster.
- This mixed state is often called partially synchronized. With further increases in K, more and more oscillators are recruited into the synchronized cluster, and $r_\infty$ tends to 1

Both of these claims are also shown to be true in the animation and in the plots of the effective velocities. In the animation, we clearly see that most of the oscillators are in a clump, but there are a few that drift away from the mean and stay out of the cluster. In the velcity plot, this is almost clearly evident. Also, if looked at carefully, it shows that the particles which are not in sync are the ones which where far away from the mean in the initial distribution"

# ╔═╡ afbd4280-8957-11eb-12a7-6762ab68086d
md"## κc vs g(ω)

The paper also analyses $\kappa_c$ as a function of $g(\omega)$. While Kuramoto and Nishikawa had previously assumed that the drifters would not have any effect on the coherence of the oscillators, they later recanted their stance and assumed that the drifters **did** in fact have an effect on the evolution of the coherence. Probably if they had the luxury of fast computers and easy computer languages, they would have seen that the drifters did in fact have an effect on $r(t)$. In the animation above, there are small but fast changes in the angle and the amplitude of the mean vector, and these are directly related to the movement of the drifters. Obviously, a simulation of only 100 oscillators is obviously not in the thermodynamic limit and the effect on $r$ is more apparent.

None the less, we can be reasonably sure that Kuramoto and Nishikawa's second attempt (§6.2 in [^strogatz]) was more accurate, and by extension so was Strogatz's own analysis. From their analysis, 

$$\kappa_c = \frac{2}{\pi g(0)}$$

In the case of the gaussian distribution, 

$$g(0) = \frac{1}{\sigma \sqrt{2\pi}}$$

Hence, 

$$\kappa_c = \frac{2\sigma\sqrt{2\pi}}{\pi} = \sqrt\frac{8}{\pi}\sigma$$
"

# ╔═╡ ff877638-8967-11eb-3c88-7b20d5b80377
md"### Computational implementation

We define a function `simulate_for_sigma` which searches by bisection for $\kappa$s which are useful for our plot. We limit the difference of $\kappa$ and difference in $r_{\infty}$ for the sake of computational time. $r_{\infty}$ is evaluated by taking a mean over the last $10\%$ timesteps.

In the ideal case, we'd average $r_\infty$ over many initial conditions to make sure that there are infact no artifacts of the initial conditions. This is especially true, near the bifurcation, where the transition times can be very long. However, we take the easier path by simply box-smoothing over the data (a tool often employed by statisticians when they get noisy data).
"

# ╔═╡ 441a57e4-889a-11eb-27b8-4bdcc4bc2826
get_rinf(arr::Array{Complex}) = mean(abs.(last(arr, 10_000)))

# ╔═╡ 1d3dcba2-8894-11eb-1429-3d90997bbc68
function simulate_for_sigma(σ::Float64, κmn::Float64, κmx::Float64, κdiff::Float64, rdiff::Float64)
	function inner_function!(d::Dict{Float64, Float64}, κmn::Float64, κmx::Float64)
		if ((κmx - κmn) < κdiff) || ((d[κmx] - d[κmn]) < rdiff)
			return
		end
		κ = (κmx + κmn)/2
		d[κ] = get_rinf(simulate(100_000, 1000, σ, κ, 100))
		inner_function!(d, κ, κmx)
		inner_function!(d, κmn, κ)
	end
	d = Dict{Float64, Float64}()
	d[κmn] = get_rinf(simulate(100_000, 1000, σ, κmn, 100))
	d[κmx] = get_rinf(simulate(100_000, 1000, σ, κmx, 100))
	inner_function!(d, κmn, κmx)
	return d
end

# ╔═╡ 8937d90e-88a0-11eb-051a-1564ca42a34a
boxsmooth(arr, b) = [mean(arr[i-b:i+b]) for i in b+1:length(arr)-b]

# ╔═╡ 42f81e66-8a78-11eb-3ee9-cbfa2e4a0cce
md"Run new simulation? $(@bind runsim CheckBox(false))"

# ╔═╡ 0778d0b8-8962-11eb-3d45-ef790f6f1701
begin
	ks = 0
	rinfs = 0
	if runsim
		d = simulate_for_sigma(3, 0., 10., 0.01, 0.01)
		ks = sort(keys(d))
		rinfs = map(ks) do k; d[k]; end
	else
		d = load("save_sigma_3_nice.jld")
		ks, rinfs = d["ks"], d["rinfs"]
	end
	plot(boxsmooth(ks, 3), boxsmooth(rinfs, 3), label="Data")
	vline!([√(8/π)*3], label="Theoretical prediction")
	vline!([1.56*3], label="\\kappa_{c} from fitting below")
	plot!(xlabel="\\kappa", ylabel="\\ r_∞", title="\\kappa vs r_∞", legend=:right)
end

# ╔═╡ 206c24ea-896b-11eb-1391-89083f16f768
md"### Stability of κc over σ

Using the powerful computers at our disposal, we can run the simulation on multiple $\kappa$s and $\sigma$s and find the $r_\infty$ for each of these. Then, we can analyse the results to find $\kappa_c$ for each $\sigma$ and fit it to a straight line and see if the slope is equal to $\sqrt{\frac{8}{\pi}}$. Like before, we smooth our data, but this time we convolve the data with a 2dimensional Gaussian blur filter.

We measure the points of phase transition by taking a derivative along $\kappa$, and picking points where the derivative is above a certain threshold(=0.09).
"

# ╔═╡ b6ce7d86-8a77-11eb-3b8e-b36f29d2839b
md"""
#### Code to generate data
```julia
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
κs, σs = collect(0.0:0.1:2.5), collect(0.0:0.1:2.5)
zs = graph(κs, σs)
save("save.jld", "zs", zs, "ks", κs, "ss", σs)
```
"""

# ╔═╡ a028f8a0-8972-11eb-05dc-3103777aa93a
begin
	kcs_data = load("save.jld")
	kcszs, kcsκs, kcsσs = kcs_data["zs"], kcs_data["ks"], kcs_data["ss"]
	kcsnzs = imfilter(kcszs, Kernel.gaussian(0.7))
	deriv = kcsnzs[2:end, :] - kcsnzs[1:end-1, :]
	j = findall(deriv .> 0.09)
	
	ys = map(j) do i; kcsκs[i[1]]; end
	xs = map(j) do i; kcsσs[i[2]]; end

	@. model(x, m) = m*x
	fit = curve_fit(model, xs, ys, [0.])
end;

# ╔═╡ d61d1470-8989-11eb-1ddb-cff31dd6d747
md"_Interactive plot below_"

# ╔═╡ 4ad38dd0-896c-11eb-2da5-335c39437209
let
	plotlyjs()
	p = plot(
		kcsκs, kcsσs, kcsnzs', 
		st=:surface, zlims=(0, 1), 
		xlabel="κ" , ylabel="σ", zlabel="r∞", 
		title="\$R_∞(\\sigma, \\kappa)\$", extra_plot_kwargs = KW(:include_mathjax => "cdn")
	)
	gr()
	p
end

# ╔═╡ e29ac97e-8989-11eb-372b-3be6756949d8
md"_Interactive plot below_"

# ╔═╡ 1f26e4dc-8973-11eb-3428-ff96adccd9d6
let
	plotlyjs()
	p = plot(
		kcsκs, kcsσs, deriv', 
		zlims=(0,0.5), st=:surface, 
		xlabel="κ" , ylabel="σ", zlabel="r∞",
		title="\$\\dfrac{dr_\\infty}{d\\kappa}(\\sigma, \\kappa)\$", extra_plot_kwargs = KW(:include_mathjax => "cdn")
	)
	gr()
	p
end

# ╔═╡ 97c0c8a4-8973-11eb-31c1-a3e993ec94b6
begin
	plot(0:0.5:2.5, model(0:0.5:2.5, fit.param), label="fit",
		yaxis="κ", xaxis="σ",
		xlims=(0,2.5), ylims=(0, 2.5))
	plot!(xs, ys, st=:scatter, label = "data",
		legend=:bottomright, title="\\kappa_{c} vs \\sigma")
end

# ╔═╡ b9a7c122-8980-11eb-1dba-a790aec65322
Markdown.parse("\$\$\\kappa_c = $(round(fit.param[1], sigdigits=3))\\cdot\\sigma\$\$")

# ╔═╡ 6b896608-898f-11eb-2c33-ff41bc3b1c87
md"### Observations

The single run and the extended analysis both show the same result that the $\kappa_c$ follows the analysis in the paper closely, however, there is a small disparity. Any finite simulation will inevitably not have the sharp nature of the bifurcation. In §4, we run the same simulations for other distributions, and see the effect of the kind of distribution on this.
"

# ╔═╡ e0a1ebd2-89a9-11eb-15b0-f5e35ecd8ffc
md"# §4 Robustness of results to distributions

As mentioned in the previous subsection, let's run the same simulation with some other distributions"

# ╔═╡ a64c11a4-89b1-11eb-2050-f5267b9c57c1
kappas = 0:0.5:10

# ╔═╡ a1d5c076-89ad-11eb-1530-cbd02ee0a46f
md"## Gaussian Distribution

We've already run this before, but let's do it again.

$$g(\omega) = \frac{1}{\sqrt{2\pi\sigma^2}}\exp\left[-\frac{\omega^2}{2\sigma^2}\right];\  g(0) = \frac{1}{\sqrt{2\pi\sigma^2}};\ \kappa_c(3) = \sqrt\frac{8}{\pi}*3 = 1.596$$
"

# ╔═╡ e88b4f7a-89af-11eb-1b9f-c9cfc4510d86
gaussian_rs = map(kappas) do kappa
	mean(abs.(last(simulate(100_000, 1000, 3., kappa, 40), 1000)))
end;

# ╔═╡ 78e8340a-89b7-11eb-219c-135ab851ccc2
begin
	plot(kappas, gaussian_rs, xlabel="\\kappa", ylabel="\\ r_∞")
	vline!([sqrt(8/pi)*3], title="Gaussian")
end

# ╔═╡ 12d46c62-89b7-11eb-1d5e-1717a0643ad1
md"## Lorenz-Cauchy distributions

$$g(\omega) = \frac{1}{\pi \sigma \left(1 + \left(\frac{\omega}{\sigma} \right)^2 \right)};\ g(0) = \frac{1}{\pi\sigma};\ \kappa_c(3) = 6$$
"

# ╔═╡ f3e77840-89b8-11eb-0d5a-b3501386bb19
begin
	cauchy_kappas = 3:0.5:12
	cauchy_rs_1 = map(cauchy_kappas) do kappa
		mean(abs.(last(simulate(50_000, init(1000, Cauchy(0,3)), kappa, 40), 1000)))
	end;
	cauchy_rs_2 = map(cauchy_kappas) do kappa
		mean(abs.(last(simulate(50_000, init(1000, Cauchy(0,3)), kappa, 40), 1000)))
	end;
	cauchy_rs_3 = map(cauchy_kappas) do kappa
		mean(abs.(last(simulate(50_000, init(1000, Cauchy(0,3)), kappa, 40), 1000)))
	end;
	cauchy_rs = (cauchy_rs_1 + cauchy_rs_2 + cauchy_rs_3)/3
end

# ╔═╡ ef8ad520-89b9-11eb-126a-db54fc2e65e9
begin
	plot(cauchy_kappas, cauchy_rs, xlabel="\\kappa", ylabel="\\ r_∞")
	vline!([6], title="Lorenz-Cauchy", ylim=(0,1))
end

# ╔═╡ 7df01b22-89c4-11eb-1422-e5aa9bd7eac6
md"## Reflected exponential decay

$$g(\omega) = \frac{1}{2a}\exp\left(-\frac{|x|}{a}\right);\ g(0) = \frac{1}{2a};\ \kappa_c(3) = \frac{4\pi}{3} = 4.189$$
"

# ╔═╡ cbdf714a-8a44-11eb-0f4d-d30ab28b7ab7
begin
	function rejectionSample(f, N)
		i, x_list = 1, zeros(N)
		while(i < N + 1)
			x=20*rand() - 10
			y=rand()
			if(f(x)>y)
				x_list[i] = x
				i += 1
			end
		end
		return x_list
	end
	red(a)= x -> exp(-abs(x)/a)/2a
end

# ╔═╡ 123f09de-8a45-11eb-13e5-fdac2199b036
red_ags = Agent.(0.0, rejectionSample(red(3), 1000));

# ╔═╡ cd5c8b66-8a49-11eb-3cf5-75946228434e
begin
	red_ags_2 = deepcopy(red_ags)
	red_rs = map(kappas) do kappa
		mean(abs.(last(simulate(100_000, red_ags_2, kappa, 40), 1000)))
	end
end;

# ╔═╡ 23872ebe-8a4b-11eb-0e7d-2d7500b874fb
begin
	plot(kappas, red_rs, xlabel="\\kappa", ylabel="\\ r_∞")
	vline!([4pi/3], title="Reflected Exponential Decay")
end

# ╔═╡ ca49a058-8a4e-11eb-1445-417745db689d
md"## Observations

Like the Gaussian, the other distributions also have a slight discrepancy from the expected value. However by a visual analysis, we can say that the theoretical results are close to what we would get in the infinite limit. For a more accurate analysis, one could fit the function as mentioned in the paper to the data, and then look at the value of $\kappa_c$.

Another observation that we can see is that distributions that are tail heavy are much more prone to variations in the partially coherent part of the graph. In order to reduce this, we average results over 3 runs with random initial conditions for the Lorenz Cauchy distribution
"

# ╔═╡ 30c5c5d4-8a57-11eb-248c-a7961f0d4590
md"# References
[^strogatz]: Steven H. Strogatz, From Kuramoto to Crawford: exploring the onset of synchronization in populations of coupled oscillators, Physica D: Nonlinear Phenomena, Volume 143, Issues 1–4, [https://doi.org/10.1016/S0167-2789(00)00094-4](https://doi.org/10.1016/S0167-2789(00)00094-4.)
[^julia]: Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah. (2017) SIAM Review, 59: 65–98. doi: [10.1137/141000671](https://dx.doi.org/10.1137/141000671)
[^plots]: Plots - powerful convenience for visualization in Julia. Thomas Breloff. [docs.juliaplots.org](docs.juliaplots.org)
[^pluto]: Pluto - Reactive Notebooks for Julia. Fons van der Plas and Mikołaj Bochenski. [https://github.com/fonsp/Pluto.jl](https://github.com/fonsp/Pluto.jl)
[^plutoui]: PlutoUI. Fons van der Plas. [https://github.com/fonsp/PlutoUI.jl](https://github.com/fonsp/PlutoUI.jl)
[^stats]: Statistics - Julia stdlib for Statistics. [https://docs.julialang.org/en/v1/stdlib/Statistics/](https://docs.julialang.org/en/v1/stdlib/Statistics/)
[^images]: JuliaImages: image processing and machine vision for Julia. [https://juliaimages.org/stable/](https://juliaimages.org/stable/)
[^lsqfit]: LsqFit.jl: Basic least-squares fitting in pure Julia under an MIT license. [https://github.com/JuliaNLSolvers/LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl)
[^distributions]: 1.Dahua Lin. JuliaStats/Distributions.jl: v0.24.15. (2021). [doi:10.5281/zenodo.4577319](https://doi.org/10.5281/zenodo.2647458)
[^jld]: JLD.jl. JuliaIO. [https://github.com/JuliaIO/JLD.jl](https://github.com/JuliaIO/JLD.jl)
"

# ╔═╡ Cell order:
# ╟─e8d31f7c-88aa-11eb-05b7-31ff25d50afb
# ╟─0d6813c6-88ae-11eb-3f39-b59ffd71ffd5
# ╟─831e6f04-8a55-11eb-35da-2f77e5dd68a5
# ╟─bd9326c2-8954-11eb-003e-b52f7107f199
# ╟─3effc012-8a56-11eb-1cd2-475e917c969b
# ╟─8d7c0eb0-84ca-11eb-0018-bdb2e389ee52
# ╟─27e115ae-88b0-11eb-1d7d-57ae925c4c9f
# ╟─bc07ce24-88b5-11eb-2052-53543e8ef9fc
# ╠═a34a7c40-84ca-11eb-3433-41615dfc1c08
# ╟─56ce942c-88b8-11eb-08fc-0970f8de6359
# ╠═c4b72810-84ca-11eb-1855-23703a18b3e4
# ╟─9a1dec4e-88e7-11eb-0736-6b6e6b24627f
# ╠═c12bde70-84ca-11eb-1fe2-634e16d69909
# ╠═9fa5bf26-89b8-11eb-0636-ed01bb9b689b
# ╟─a7d91872-88e7-11eb-2799-3502fc01e720
# ╠═cb3e73a0-84ca-11eb-3c68-8fc3fc5cf239
# ╠═d00b64fe-84ca-11eb-1e68-f73bbca9e373
# ╟─393ed800-88e9-11eb-2160-edef4c90ee74
# ╟─86408954-88f9-11eb-2ee7-85532daf659d
# ╠═8c9b5c8e-84d3-11eb-0f9d-47c844a46a68
# ╟─c618ba40-84d2-11eb-3990-93a786c8f8a9
# ╟─cbb28294-88f9-11eb-1b77-fdacf61b8c9b
# ╟─f564ad60-88ea-11eb-31c5-a9f36ba1c93b
# ╠═bef4689a-88ec-11eb-2297-3bfe711d504b
# ╟─2a5fcef8-88ed-11eb-36e9-a743d38881f6
# ╟─586e0782-88f3-11eb-2f30-03211c7a70a2
# ╟─f2a5c968-88fb-11eb-3a39-2b95d2a8d978
# ╟─222d4b00-88f4-11eb-2f5b-bd320a59b2ab
# ╟─9537f2f8-88f4-11eb-091d-b5827649cb40
# ╟─def3ccb4-8949-11eb-1ab1-a3403dc946ad
# ╟─d2d048a8-8956-11eb-0536-972068f6862f
# ╟─963c13ba-8957-11eb-3d45-f1d14f7b1c87
# ╟─afbd4280-8957-11eb-12a7-6762ab68086d
# ╟─ff877638-8967-11eb-3c88-7b20d5b80377
# ╠═441a57e4-889a-11eb-27b8-4bdcc4bc2826
# ╠═1d3dcba2-8894-11eb-1429-3d90997bbc68
# ╠═8937d90e-88a0-11eb-051a-1564ca42a34a
# ╟─42f81e66-8a78-11eb-3ee9-cbfa2e4a0cce
# ╠═0778d0b8-8962-11eb-3d45-ef790f6f1701
# ╟─206c24ea-896b-11eb-1391-89083f16f768
# ╟─b6ce7d86-8a77-11eb-3b8e-b36f29d2839b
# ╠═a028f8a0-8972-11eb-05dc-3103777aa93a
# ╟─d61d1470-8989-11eb-1ddb-cff31dd6d747
# ╟─4ad38dd0-896c-11eb-2da5-335c39437209
# ╟─e29ac97e-8989-11eb-372b-3be6756949d8
# ╟─1f26e4dc-8973-11eb-3428-ff96adccd9d6
# ╟─97c0c8a4-8973-11eb-31c1-a3e993ec94b6
# ╟─b9a7c122-8980-11eb-1dba-a790aec65322
# ╟─6b896608-898f-11eb-2c33-ff41bc3b1c87
# ╟─e0a1ebd2-89a9-11eb-15b0-f5e35ecd8ffc
# ╠═a64c11a4-89b1-11eb-2050-f5267b9c57c1
# ╟─a1d5c076-89ad-11eb-1530-cbd02ee0a46f
# ╠═e88b4f7a-89af-11eb-1b9f-c9cfc4510d86
# ╟─78e8340a-89b7-11eb-219c-135ab851ccc2
# ╟─12d46c62-89b7-11eb-1d5e-1717a0643ad1
# ╠═f3e77840-89b8-11eb-0d5a-b3501386bb19
# ╟─ef8ad520-89b9-11eb-126a-db54fc2e65e9
# ╟─7df01b22-89c4-11eb-1422-e5aa9bd7eac6
# ╟─cbdf714a-8a44-11eb-0f4d-d30ab28b7ab7
# ╠═123f09de-8a45-11eb-13e5-fdac2199b036
# ╠═cd5c8b66-8a49-11eb-3cf5-75946228434e
# ╟─23872ebe-8a4b-11eb-0e7d-2d7500b874fb
# ╟─ca49a058-8a4e-11eb-1445-417745db689d
# ╟─30c5c5d4-8a57-11eb-248c-a7961f0d4590
