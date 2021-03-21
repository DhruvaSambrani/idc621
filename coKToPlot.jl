using JLD
using Plots
using Images
using LsqFit

plotlyjs()

d = load("save.jld")
zs, κs, σs = d["zs"], d["ks"], d["ss"]
nzs = imfilter(zs, Kernel.gaussian(0.7))
plot(κs, σs, nzs', st=:surface, zlims=(0, 1), xlabel="κ" , ylabel="σ", zlabel="r_inf", title="R∞(σ, κ)")
display(current())

k = nzs[2:end, :] - nzs[1:end-1, :]

plot(κs, σs, k', zlims=(0,0.5), st=:surface, xlabel="κ" , ylabel="σ", zlabel="r_inf", title="R∞(σ, κ)")
display(current())

j = findall(k .> 0.09)

ys = map(j) do i; κs[i[1]]; end
xs = map(j) do i; σs[i[2]]; end

@. model(x, m) = m*x
fit = curve_fit(model, xs, ys, [0.])

plot(0:0.5:2.5, model(0:0.5:2.5, fit.param), yaxis="κ", xaxis="σ", xlims=(0,2.5), ylims=(0, 2.5))
plot!(xs, ys, st=:scatter, legend=:bottomright)

display(current())
