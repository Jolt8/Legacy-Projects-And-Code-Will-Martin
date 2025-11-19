using ModelingToolkit
using DifferentialEquations
using Plots
using GLMakie

import ModelingToolkit: t_nounits as t, D_nounits as D

using Ferrite

#this was taken from the MTK standard library 
@connector HeatPort begin
    @parameters begin
        T_guess = 273.15 + 20
        Q_flow_guess = 0.0
    end

    @variables begin
        T(t), [guess = T_guess]
        Q_flow(t), [guess = Q_flow_guess, connect = Flow]
    end
end

# A simple 0D heat capacitor, representing one node.
@mtkmodel HeatCapacitor begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        C = 1.0, [description = "Heat capacity of element"]
    end
    @variables begin
        T(t) = 293.15, [description = "Temperature of element [K]"]
        der_T(t) = 0.0, [description = "Time derivative of temperature [K/s]"]
    end

    @equations begin
        T ~ port.T
        D(T) ~ port.Q_flow / C
    end
end

@mtkmodel Element1D begin
    @components begin
        port_a = HeatPort()
        port_b = HeatPort()
    end
    @variables begin
        dT(t), [guess = 0.0]
        Q_flow(t), [guess = 0.0]
    end
    @equations begin
        dT ~ port_a.T - port_b.T
        port_a.Q_flow ~ Q_flow
        port_a.Q_flow + port_b.Q_flow ~ 0
    end
end

# A 0D thermal conductor, representing the connection between nodes.
@mtkmodel ThermalConductor begin
    @extend Q_flow, dT = element1d = Element1D()
    @parameters begin
        G = 1.0, [description = "Thermal conductance [W/K]"]
    end
    @equations begin
        Q_flow ~ G * dT
    end
end

# A constant temperature boundary condition
@mtkmodel SomeFixedTemperature begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        T_fixed = 303.15, [description = "Fixed temperature [K]"] # 30 °C
    end
    @equations begin
        port.T ~ T_fixed
    end
end

# A constant heat flow boundary condition
@mtkmodel SomeHeatSource begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        Q_flow_fixed = 10.0, [description = "Fixed heat flow [W]"]
    end
    @equations begin
        port.Q_flow ~ Q_flow_fixed
    end
end

#START OF FERRITE SECTION
grid = generate_grid(Quadrilateral, (1, 1))

addnodeset!(grid, "left", (x) -> x[1] ≈ -1)
#getnodeset(grid, "left")
addnodeset!(grid, "right", (x) -> x[1] ≈ 1)

#collect(getnodeset(grid, "left"))
#collect(getnodeset(grid, "right"))

top = ExclusiveTopology(grid)

vertex_neighbors = Ferrite.vertex_star_stencils(top, grid)

unique_connections = Set{Tuple{Int, Int}}()

for i in eachindex(vertex_neighbors)
    node_at_i = collect(vertex_neighbors[i][1].idx)[2]
    for j in 2:length(vertex_neighbors[i])
        neighbor_at_i = collect(vertex_neighbors[i][j].idx)[2]
        #println(i, ", ", j, ", ", vertex_neighbors[i][j])
        #println(node_at_i, ", ", neighbor_at_i)
        #println((min(node_at_i, neighbor_at_i), max(node_at_i, neighbor_at_i)))
        push!(unique_connections, (min(node_at_i, neighbor_at_i), max(node_at_i, neighbor_at_i)))
    end
end

#println(unique_connections)

unique_connections = collect(unique_connections)
#END OF FERRITE SECION

n_nodes = getnnodes(grid)

nodes = [HeatCapacitor(name=Symbol("node_", i), C=1/n_nodes) for i in 1:n_nodes]

conds = [ThermalConductor(name=Symbol("h_cond_", i)) for i in 1:length(unique_connections)]

#Add boundary condition components
#thought this would be hard but getnodeset(grid, "some_boundary_condition_name") does exactly what I need
n_left_bcs = length(collect(getnodeset(grid, "left")))
left_bcs = [SomeHeatSource(name=Symbol("heat_source_", i)) for i in 1:n_left_bcs] 
n_right_bcs = length(collect(getnodeset(grid, "right")))
right_bcs = [SomeFixedTemperature(T_fixed=293.15, name=Symbol("fixed_temp_", i)) for i in 1:n_right_bcs]

connections = Equation[]

for (i, (i_in, i_out)) in enumerate(unique_connections)
    push!(connections, connect(nodes[i_in].port, conds[i].port_a))
    push!(connections, connect(conds[i].port_b, nodes[i_out].port))
end

left_node_indicies = collect(getnodeset(grid, "left"))
for (i, node_idx) in enumerate(left_node_indicies)
    push!(connections, connect(left_bcs[i].port, nodes[node_idx].port))
end

right_node_indicies = collect(getnodeset(grid, "right"))
for (i, node_idx) in enumerate(right_node_indicies)
    push!(connections, connect(right_bcs[i].port, nodes[node_idx].port))
end

# START OF SIM SETUP
eqs_total = System[]
all_systems = vcat(
    vec(nodes),
    vec(conds),
    vec(left_bcs),
    vec(right_bcs)
)
all_systems

@named rod = ODESystem(connections, t, systems=all_systems)

@time begin 
    sys = structural_simplify(rod)  
    
    tspan = (0.0, 10.0) #changing the tspan does not cause a second run long run time

    prob = ODEProblem(sys, [], tspan)

    sol = solve(prob)
end

#if you want to time a file: @time include("C://Users//wille//Desktop//Legacy-Projects-And-Code-Will-Martin//Julia Projects//Testing FEM with mtkcomponents.jl")

#these times are all after the program has compiled 
#the first time is when its run with a new grid shape/size the second time is after rerunning with the same grid shape/size
#this is tested on a i7-7700HQ

#How much does first time compilation take on a 1x1 grid vs a 10x10 grid
#if the compilation time (not the total time to run the file) is less on the 1x1 grid, would it be a good idea to initialize it with the simplest grid possible 
#then run a more complex grid?

#1x1 takes 10.801329 seconds (23.36 M allocations: 1.209 GiB, 3.84% gc time, 98.04% compilation time: <1% of which was recompilation)
#1x1 takes 0.225861 seconds (806.67 k allocations: 27.395 MiB, 0.01% compilation time: 100% of which was recompilation)

#just as an example 1x2 from cold start takes 64.063055 seconds (112.72 M allocations: 5.789 GiB, 4.54% gc time, 99.43% compilation time: 46% of which was recompilation)
#1x2 takes 14.131481 seconds (24.14 M allocations: 1.236 GiB, 18.59% gc time, 97.39% compilation time: <1% of which was recompilation)
#1x2 takes 0.462619 seconds (1.41 M allocations: 50.913 MiB, 0.01% compilation time: 100% of which was recompilation)

#5x5 takes 16.7 seconds at first
#5x5 takes 1.21 seconds at second

#10x10 takes 25.2 seconds at first
#1x1 takes 2.83 seconds at second

#15x15 takes 53.8 seconds at first
#15x15 takes 9.05 seconds at second

#I'm starting to realize that this isn't very scalable, god damn it

#20x20 takes 134.05 seconds at first
#20x20 takes 17.42 seconds at second

#30x30 takes 448.80 seconds at first
#30x30 takes 104.98 seconds at second (149.42 M allocations: 13.025 GiB, 6.27% gc time, 0.00% compilation time: 100% of which was recompilation)

#30x60 takes 1220.276316 seconds (557.49 M allocations: 48.750 GiB, 2.27% gc time, 71.29% compilation time: <1% of which was recompilation)
#30x60 takes 305.805608 seconds (457.29 M allocations: 44.694 GiB, 3.92% gc time, 0.00% compilation time: 94% of which was recompilation)


#I wonder what takes more time a 1x16, 2x8, 4x4 grid

#1x16 takes 19.260376 seconds (28.80 M allocations: 1.448 GiB, 12.91% gc time, 93.78% compilation time: <1% of which was recompilation)
#1x16 takes 1.363001 seconds (4.30 M allocations: 199.238 MiB, 12.31% gc time, 0.01% compilation time: 100% of which was recompilation)

#2x8 takes 16.647645 seconds (26.93 M allocations: 1.345 GiB, 2.89% gc time, 94.49% compilation time: <1% of which was recompilation)
#2x8 takes 1.337417 seconds (3.17 M allocations: 124.089 MiB, 0.01% compilation time: 100% of which was recompilation)

#4x4 takes 20.704495 seconds (26.23 M allocations: 1.313 GiB, 3.16% gc time, 95.33% compilation time: <1% of which was recompilation)
#4x4 takes 1.042183 seconds (2.70 M allocations: 97.871 MiB, 0.01% compilation time: 100% of which was recompilation)

#the good news is that even if a grid that was previously compiled is used again after being replaced it will not have to 
#recompile and you'll get the 100% recompilation time
#ex. setting it to (1, 16) running it, changing to (2, 8), running (1, 16) again preserves compilation

#also, usually the time it takes to recompile (not simulation) goes down slightly every time you run it

#also, changing a parameter like C of each heat capacitor does not require another simulation (only recompilation) which is actually REALLY good

#after looking at this in sheets the first time to run seems to be T(n log n) and the second time to run seems like halfway between T(n log n) and T(n^2)
#Here's the link to the sheet: https://docs.google.com/spreadsheets/d/1nOqnXJgxIux2OHeiS7nu_4JaD5Po5ocJ0T5oBaHbS3A/edit?usp=sharing


