# %%
include("CellularBase.jl")
using .CellularBase
println("using CB")
@time using Plots
println("using Plots")
@time using Images
println("using Images")
# %%
mutable struct LifeGrid <: AbstractGrid{Int}
    state::Array{Int}
    otherstate::Array{Int}
    bc::BoundaryCondition
    neighborhood::AbstractNeighborhood
    τi::Int64
    τr::Int64
    max::Float64
    min::Float64
    function LifeGrid(
            size::Tuple{Int,Int}, τi::Int, τr::Int;
            max::Float64, min::Float64, init=nothing)
        z = zeros(size)
        if isnothing(init) z[(size .÷ 2)...] = 1
        else z[rand(CartesianIndices(z), init)] .= 1
        end
        new(z, copy(z), FixedMax, VonNeumannNeighborhood(3, 2), τi, τr, max, min)
    end
end

CellularBase.possible_states(l::LifeGrid) = 0:(l.τi + l.τr)

CellularBase.newstate(grid::LifeGrid) = grid.otherstate;

function CellularBase.state!(grid::LifeGrid, newstate)
    grid.otherstate = grid.state
    grid.state = newstate
end
# %%
infect(fraction::Float64, g::LifeGrid) = rand() < max(((g.max - g.min) * fraction + g.min), 0)
susceptible(i::Int, g::LifeGrid) = i == 0
infectious(i::Int, g::LifeGrid) = i ∈ 1:g.τi
refactory(i::Int, g::LifeGrid) = i ∈ g.τi:(g.τr + g.τi)
recovered(i::Int, g::LifeGrid) = i == (g.τr + g.τi)
function (g::LifeGrid)(neighbors::Array{Int}; kwargs...)
	current = neighbors[1]
    if susceptible(current, g) 
        fraction = count(i -> infectious(i, g), neighbors) / length(neighbors)
        if infect(fraction, g)
            return 1
        else
            return 0
        end
    else
        return (current + 1) % (g.τi + g.τr)
    end
end
function color(i, grid)
    if susceptible(i, grid)
        RGB(0.9, 0.9, 0.9)
    elseif infectious(i, grid)
        RGB(1., 0., 0.)
    else
        RGB(0., 0., 0.)
    end
end;
# %%
function lockdown(grid, step; step_dict, radii)
    push!(radii, radius(grid.neighborhood))
    if step in keys(step_dict)
        grid.neighborhood = VonNeumannNeighborhood(step_dict[step], 2); 
    end
end
step_dict = Dict(50 => 2, 100 => 1, 200 => 2, 300 => 3)
radii = []
grid = LifeGrid((200, 200), 4, 9, min=0., max=0.5, init=5)
k = simulate!(grid, 400; postrunhook=lockdown, step_dict=step_dict, radii=radii);
# %%
begin
    gif = @gif for frame ∈ k
        plot(frame .|> p -> color(p, grid))
    end every 2
    display(gif)
end
# %%
@time begin
    p = plot(1:length(k), k .|> frame -> count(frame .|> p -> susceptible(p, grid)), color=:green, label="Sus")
    plot!(1:length(k), k .|> frame -> count(frame .|> p -> infectious(p, grid)), color=:red, label="Inf")
    plot!(1:length(k), k .|> frame -> count(frame .|> p -> refactory(p, grid)), color=:blue, label="Ref")
    q = plot(radii, label="Neighbourhood Radius", legend=:bottomright)
    plot(p, q, size=(600, 500), layout=(2, 1))
end
# %%
