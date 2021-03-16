using JLD
using Plots
using Images
using LsqFit

plotlyjs()

d = load("save.jld")
zs, κs, σs = d["zs"], d["ks"], d["ss"]

plot(σs, κs, zs, st=:surface, zlims=(0, 1), xlabel="Sigma", ylabel="Kappa", zlabel="r_inf", title="R∞(σ, κ)")
display(current())

nzs = imfilter(zs, Kernel.gaussian(0.7))

plot(σs, κs, nzs, st=:surface, zlims=(0, 1), xlabel="Sigma", ylabel="Kappa", zlabel="r_inf", title="R∞(σ, κ)")
display(current())

k = nzs[:, 1:end-1] - nzs[:, 2:end]

plot(k, st=:wireframe)
display(current())

j = findall(k .> 0.2)

plot(k .> 0.2, st=:wireframe)
display(current())

xs = map(j) do i; κs[i[1]]; end
ys = map(j) do i; σs[i[2]]; end

@. model(x, (m,c)) = m*x+c

fit = curve_fit(model, xs, ys, [0.,0.])

plot(0:0.5:2, model(0:0.5:2, fit.param), xaxis="κ", yaxis="Sigma")
plot!(xs, ys, st=:scatter, legend=:bottomright)

display(current())
