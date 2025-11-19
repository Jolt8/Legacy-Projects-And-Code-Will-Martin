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
        k = 1.0, [description = "thermal conductivity [W/(m*k)]"]
        A = 1.0, [description = "facet area [m^2]"]
        d = 1.0, [description = "distance between centroids [m]"]
        C = k * (A / d), [description = "Heat capacity of element [J/(g*K)]"]
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

grid = generate_grid(Hexahedron, (2, 2, 2), left, right)

addcellset!(grid, "left", (x) -> x[1] ≈ 0)
addcellset!(grid, "right", (x) -> x[1] ≈ 1)

function generate_FVM_connectivity(grid, poly_interp, cell_quadrature_rule, facet_quadrature_rule)
    cell_qr = cell_quadrature_rule
    cell_values = CellValues(cell_qr, poly_interp)
    facet_qr = facet_quadrature_rule
    facet_values = FacetValues(facet_qr, poly_interp)

    cell_dimensions = []

    for cell in CellIterator(grid)
        #Get cell volume
        reinit!(cell_values, cell)
        vol_i = 0.0 
        for qp in 1:getnquadpoints(cell_values)
            d_vol = getdetJdV(cell_values, qp) 
            vol_i += d_vol
        end
        
        #Get face areas
        area = []
        cell_coords = getcoordinates(cell)
        
        for facet_index in 1:nfacets(cell)
            reinit!(facet_values, cell_coords, facet_index)
            
            area_i = 0.0
            for qp in 1:getnquadpoints(facet_values)
                area_i += getdetJdV(facet_values, qp)
            end
            
            push!(area, area_i)
        end

        #Get centroid
        x_avg = sum([c[1] for c in cell_coords]) / length(cell_coords)
        y_avg = sum([c[2] for c in cell_coords]) / length(cell_coords)
        z_avg = sum([c[3] for c in cell_coords]) / length(cell_coords)

        centroid_coords = [x_avg, y_avg, z_avg]

        push!(cell_dimensions, [vol_i, [area], centroid_coords]) #formatting it this way because it's kindof 3D, 2D, 1D although a different format might be better
    end
    return cell_dimensions
end

poly_interp = Lagrange{RefHexahedron, 1}() #1 = linear elements 2 = quadratic/curved edges

cell_qr = QuadratureRule{RefHexahedron}(2) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation

facet_qr = FacetQuadratureRule{RefHexahedron}(2) 

cell_dimensions = generate_FVM_connectivity(grid, poly_interp, cell_qr, facet_qr)
#END OF FERRITE SECION

n_cells = getncells(grid)

cells = [HeatCapacitor(name=Symbol("node_", i), C=cell_dimensions[i][1]*1.7*1000000) for i in 1:n_cells] #get volume #using a test value of 1.7 J/(g*K) and 1000000g/m^3

conds = [ThermalConductor(name=Symbol("h_cond_", i)) for facet_idx in 1:nfacets(getcells(grid)[i]) for cell_idx in 1:n_cells] # need a way to multiply U by the area between the cells
#also, won't this triple up thermal conductors because each should really define 3 connections to other nodes
#I guess we need a similar unique_connections to the FEM example

#Add boundary condition components
#thought this would be hard but getnodeset(grid, "some_boundary_condition_name") does exactly what I need
n_left_bcs = length(collect(getcellset(grid, "left")))
left_bcs = [SomeHeatSource(name=Symbol("heat_source_", i)) for i in 1:n_left_bcs] 
n_right_bcs = length(collect(getcellset(grid, "right")))
right_bcs = [SomeFixedTemperature(T_fixed=293.15, name=Symbol("fixed_temp_", i)) for i in 1:n_right_bcs]

connections = Equation[]

#START of small ferrite section
function FVM_unique_connections(grid)
    top = ExclusiveTopology(grid)
    unique_connections = Set() #Set{Tuple{Tuple{Number, Number}, Tuple{Number, Number}}}() #don't know if not defning the set type decreases performance
    for cell in CellIterator(grid)
        i = cellid(cell) 
        for j in 1:nfacets(cell)
            neighbors = top.face_face_neighbor[i, j]
            if !isempty(neighbors)
                push!(unique_connections, minmax((i, j), neighbors[1].idx))
            end
        end
    end
    return unique_connections
end

unique_connections = FVM_unique_connections(grid)
#END of small ferrite section
#unique_connections is formatted [[[cell_idx, ]...]

for (i, ((cell_idx, facet_idx), (neighbor_idx, neighbor_facet_idx))) in enumerate(unique_connections)
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

sys = structural_simplify(rod)  

tspan = (0.0, 10.0)

prob = ODEProblem(sys, [], tspan)

sol = solve(prob)
