using ModelingToolkit
using DifferentialEquations
#using Plots
#using GLMakie

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
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

grid_dimensions = (10, 10, 1)

grid = generate_grid(Hexahedron, grid_dimensions, left, right)

length_to_node_ratio = right[1] / collect(grid_dimensions)[1]

addnodeset!(grid, "left", x -> x[1] <= left[1] + length_to_node_ratio)
getnodeset(grid, "left")
addnodeset!(grid, "right", (x) -> x[1] >= (right[1] - 0.0000001) - (length_to_node_ratio)) #1 doesn't work
getnodeset(grid, "right")

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

unique_connections

for (i, (i_in, i_out)) in enumerate(unique_connections)
    push!(connections, connect(nodes[i_in].port, conds[i].port_a))
    push!(connections, connect(conds[i].port_b, nodes[i_out].port))
end

left_node_indicies = collect(getnodeset(grid, "left"))
for (i, node_idx) in enumerate(left_node_indicies)
    push!(connections, connect(left_bcs[i].port, nodes[node_idx].port))
end
"""
right_node_indicies = collect(getnodeset(grid, "right"))
for (i, node_idx) in enumerate(right_node_indicies)
    push!(connections, connect(right_bcs[i].port, nodes[node_idx].port))
end
"""
# START OF SIM SETUP
eqs_total = System[]
all_systems = vcat(
    vec(nodes),
    vec(conds),
    vec(left_bcs),
    #vec(right_bcs)
) #I'm starting to realize that the right_bcs wasn't necessary as it makes the system fully determined... oops.
all_systems

@named rod = ODESystem(connections, t, systems=all_systems)

@time begin 
    sys = structural_simplify(rod)  
    
    tspan = (0.0, 10.0) #changing the tspan does not cause a second run long run time

    prob = ODEProblem(sys, [], tspan)

    sol = solve(prob)
end

#look at the non 3D version of this file for timing information for 2D grids

#3x3x3 takes 27.438133 seconds (39.10 M allocations: 1.883 GiB, 2.43% gc time, 89.26% compilation time: <1% of which was recompilation)
#3x3x3 takes 3.516172 seconds (10.68 M allocations: 495.219 MiB, 6.29% gc time, 0.17% compilation time: 100% of which was recompilation)

#2x2x2 with right bcs takes 18.595824 seconds (30.95 M allocations: 1.522 GiB, 9.54% gc time, 83.59% compilation time: <1% of which was recompilation)
#2x2x2 with right bcs takes 1.602803 seconds (6.17 M allocations: 277.451 MiB, 6.36% gc time, 0.01% compilation time: 100% of which was recompilation)

#2x2x2 without right bcs takes 11.350268 seconds (14.44 M allocations: 662.630 MiB, 2.88% gc time, 81.70% compilation time: <1% of which was recompilation)
#2x2x2 without right bcs takes 2.025058 seconds (7.23 M allocations: 335.966 MiB, 9.59% gc time, 0.01% compilation time: 100% of which was recompilation)
#why does the one without right bcs take longer than the one with right bcs on their second runs
#this might mean that second time running scales differently than I thought

#4x4x4 without right bcs takes 49.184115 seconds (48.57 M allocations: 2.477 GiB, 2.73% gc time, 88.69% compilation time: <1% of which was recompilation)
#4x4x4 without right bcs takes 6.267527 seconds (17.27 M allocations: 1006.160 MiB, 28.73% gc time, 0.00% compilation time: 100% of which was recompilation)