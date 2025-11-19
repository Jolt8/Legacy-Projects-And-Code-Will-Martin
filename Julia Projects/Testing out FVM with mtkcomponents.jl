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
        #k = 1.0, [description = "thermal conductivity [W/(m*k)]"]
        #A = 1.0, [description = "facet area [m^2]"]
        #d = 1.0, [description = "distance between centroids [m]"]
        C = 1.0 
        #k * (A / d), [description = "Heat capacity of element [J/(g*K)]"]
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
        #Q_flow ~ G * dT
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

grid = generate_grid(Hexahedron, (4, 1, 4), left, right)

addcellset!(grid, "left", x -> x[1] <= 0.25)
getcellset(grid, "left")
addcellset!(grid, "right", (x) -> x[1] >= 0.75)
getcellset(grid, "right")

struct CellGeometry
    volume::Float64
    face_areas::Vector{Float64}
    centroid_coords::Vector{Float64}
    unique_connections::Vector{Any}
end

function FVM_cell_geometries(grid, poly_interp, cell_quadrature_rule, facet_quadrature_rule)
    cell_qr = cell_quadrature_rule
    cell_values = CellValues(cell_qr, poly_interp)
    facet_qr = facet_quadrature_rule
    facet_values = FacetValues(facet_qr, poly_interp)

    cell_geometries = Vector{CellGeometry}()

    for cell in CellIterator(grid)
        #Get cell volume
        Ferrite.reinit!(cell_values, cell)
        vol = 0.0 
        for qp in 1:getnquadpoints(cell_values)
            d_vol = getdetJdV(cell_values, qp) 
            vol += d_vol
        end
        
        #Get face areas
        areas = []
        cell_coords = getcoordinates(cell)
        
        for facet_index in 1:nfacets(cell)
            Ferrite.reinit!(facet_values, cell_coords, facet_index)
            
            area_i = 0.0
            for qp in 1:getnquadpoints(facet_values)
                area_i += getdetJdV(facet_values, qp)
            end
            
            push!(areas, area_i)
        end

        #Get centroid
        centroid_vec = sum(cell_coords) / length(cell_coords)

        centroid_coords = [centroid_vec[1], centroid_vec[2], centroid_vec[3]]

        push!(cell_geometries, CellGeometry(vol, areas, centroid_coords, [])) #formatting it this way because it's kindof 3D, 2D, 1D although a different format might be better
    end
    
    top = ExclusiveTopology(grid)
    unique_connections = Set{Vector{Vector{Int}}}() #Set{Tuple{Tuple{Number, Number}, Tuple{Number, Number}}}() 

    for cell in CellIterator(grid)
        i = cellid(cell) 
        for j in 1:nfacets(cell)
            neighbors = top.face_face_neighbor[i, j]
            if !isempty(neighbors)
                push!(unique_connections, collect(minmax([i, j], collect(neighbors[1].idx))))
            end
        end
    end

    println("HI", unique_connections)

    unique_connections = collect(unique_connections)

    for i in 1:length(unique_connections) #should probably find a way to connect this with the loop above but I don't know how without a set to get rid of duplicates
        cell_idx = unique_connections[i][1][1]
        println(cell_idx)
        push!(cell_geometries[cell_idx].unique_connections, unique_connections[i])
    end
    return cell_geometries
end

poly_interp = Lagrange{RefHexahedron, 1}() #1 = linear elements 2 = quadratic/curved edges

cell_qr = QuadratureRule{RefHexahedron}(2) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation

facet_qr = FacetQuadratureRule{RefHexahedron}(2) 

cell_geometries = FVM_cell_geometries(grid, poly_interp, cell_qr, facet_qr)

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
#END OF FERRITE SECION

k_thermal = 200
rho = 2700
cp = 900
vol_heat_cap = rho * cp

n_cells = getncells(grid)

cells = [HeatCapacitor(name=Symbol("node_", i), C=cell_geometries[i].volume*vol_heat_cap) for i in 1:n_cells] 

conds = []

for i in 1:length(cell_geometries)
    for j in 1:length(cell_geometries[i].unique_connections)
        face_area = cell_geometries[i].face_areas[cell_geometries[i].unique_connections[j][1][2]]
        neighbor_cell_idx = cell_geometries[i].unique_connections[j][2][2]
        x_dist = (cell_geometries[i].centroid_coords[1] + cell_geometries[neighbor_cell_idx].centroid_coords[1])
        y_dist = (cell_geometries[i].centroid_coords[2] + cell_geometries[neighbor_cell_idx].centroid_coords[2])
        z_dist = (cell_geometries[i].centroid_coords[3] + cell_geometries[neighbor_cell_idx].centroid_coords[3])
        dist_between_cells = x_dist + y_dist + z_dist #not accurate but good enough for now
        G_calc = k_thermal * face_area / dist_between_cells
        push!(conds, ThermalConductor(name=Symbol("h_cond_", i, "_", j), G=cell_geometries[i].face_areas[j]), )
    end
end

#conds = [ThermalConductor(name=Symbol("h_cond_")), G=cell.geometries[i] for cell_geom in cell_geometries for j in 1:length(cell_geom.unique_connections)] # need a way to multiply U by the area between the cells

#Add boundary condition components
#thought this would be hard but getnodeset(grid, "some_boundary_condition_name") does exactly what I need
n_left_bcs = length(collect(getcellset(grid, "left")))
left_bcs = [SomeHeatSource(name=Symbol("heat_source_", i), Q_flow_fixed=5000) for i in 1:n_left_bcs] 
n_right_bcs = length(collect(getcellset(grid, "right")))
right_bcs = [SomeFixedTemperature(name=Symbol("fixed_temp_", i), T_fixed=293.15) for i in 1:n_right_bcs]

connections = Equation[]

for i in 1:length(cell_geometries)
    for j in 1:length(cell_geometries)
        push!(connections, connect(cells[cell_idx].port, conds[facet_idx].port_a))
        push!(connections, connect(conds[facet_idx].port_b, cells[neighbor_idx].port))
    end
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
