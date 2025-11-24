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


grid_dimensions = (2, 2, 2)

#grid_dimensions = (3, 3, 3)

#grid_dimensions = (4, 4, 4)

#grid_dimensions = (5, 5, 5)


grid = generate_grid(Hexahedron, grid_dimensions, left, right)

length_to_node_ratio = right[1] / collect(grid_dimensions)[1]

addcellset!(grid, "left", x -> x[1] <= left[1] + length_to_node_ratio)
getcellset(grid, "left")
addcellset!(grid, "right", (x) -> x[1] >= (right[1] - 0.0000001) - (length_to_node_ratio)) #1 doesn't work
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

    unique_connections = collect(unique_connections)

    for i in 1:length(unique_connections) #should probably find a way to connect this with the loop above but I don't know how without a set to get rid of duplicates
        cell_idx = unique_connections[i][1][1]
        push!(cell_geometries[cell_idx].unique_connections, unique_connections[i])
    end
    return cell_geometries
end

poly_interp = Lagrange{RefHexahedron, 1}() #1 = linear elements 2 = quadratic/curved edges

cell_qr = QuadratureRule{RefHexahedron}(2) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation

facet_qr = FacetQuadratureRule{RefHexahedron}(2) 

cell_geometries = FVM_cell_geometries(grid, poly_interp, cell_qr, facet_qr)

#END OF FERRITE SECION

k_thermal = 200
rho = 2700
cp = 900
vol_heat_cap = rho * cp

n_cells = getncells(grid)

cells = [HeatCapacitor(name=Symbol("node_", i), C=cell_geometries[i].volume*vol_heat_cap) for i in 1:n_cells] 

conds = []

for i in 1:length(cell_geometries)
    push!(conds, [])
    for j in 1:length(cell_geometries[i].unique_connections)
        face_area = cell_geometries[i].face_areas[cell_geometries[i].unique_connections[j][1][2]]
        neighbor_cell_idx = cell_geometries[i].unique_connections[j][2][2]
        x_dist = (cell_geometries[i].centroid_coords[1] + cell_geometries[neighbor_cell_idx].centroid_coords[1])
        y_dist = (cell_geometries[i].centroid_coords[2] + cell_geometries[neighbor_cell_idx].centroid_coords[2])
        z_dist = (cell_geometries[i].centroid_coords[3] + cell_geometries[neighbor_cell_idx].centroid_coords[3])
        dist_between_cells = x_dist + y_dist + z_dist #not accurate but good enough for now

        dist_between_cells = norm(collect(cell_geometries[i].centroid_coords) - collect(cell_geometries[neighbor_cell_idx].centroid_coords))
        G_calc = k_thermal * face_area / dist_between_cells
        push!(conds[i], ThermalConductor(name=Symbol("h_cond_", i, "_", j), G=cell_geometries[i].face_areas[j]))
    end
end

#Add boundary condition components
left_node_indicies = collect(getcellset(grid, "left"))
n_left_bcs = length(left_node_indicies)
left_bcs = [SomeHeatSource(name=Symbol("heat_source_", i), Q_flow_fixed=5000) for i in 1:n_left_bcs] 
right_node_indicies = collect(getcellset(grid, "right"))
n_right_bcs = length(right_node_indicies)
right_bcs = [SomeFixedTemperature(name=Symbol("fixed_temp_", i), T_fixed=293.15) for i in 1:n_right_bcs]

connections = Equation[]
#cell_geometries[1].unique_connections

for i in 1:length(cell_geometries)
    for j in 1:length(cell_geometries[i].unique_connections)
        push!(connections, connect(cells[i].port, conds[i][j].port_a))
        push!(connections, connect(conds[i][j].port_b, cells[cell_geometries[i].unique_connections[j][2][1]].port))
    end
end

for (i, cell_idx) in enumerate(left_node_indicies)
    push!(connections, connect(left_bcs[i].port, cells[cell_idx].port))
end
"""
for (i, cell_idx) in enumerate(right_node_indicies)
    push!(connections, connect(right_bcs[i].port, cells[cell_idx].port))
end
"""
new_conds = []
for i in 1:length(cell_geometries)
    for j in 1:length(cell_geometries[i].unique_connections)
        push!(new_conds, conds[i][j])
    end
end

# START OF SIM SETUP
eqs_total = System[]
all_systems = vcat(
    vec(cells),
    vec(new_conds),
    vec(left_bcs),
    #vec(right_bcs)
)

#@timed begin
@named rod = ODESystem(connections, t, systems=all_systems)

println("structural_simplify time")
@time sys = structural_simplify(rod) #oh wow, structural simplify only takes 7.6 seconds on a 5x5x5 grid
#not what's better, structural simplify vs mtkcompile vs mtkbuild

tspan = (0.0, 10.0)

println("prob = ODEProblem() time")
@time prob = ODEProblem(sys, [], tspan) #it seems like this take a majority of the time, a 5x5x5 grid takes 129.5 seconds

#=
Apparently if we need to use the 2nd run time speed again we can do this
    #using Serialization
    #serialize("my_compiled_problem.jls", prob)

Then in another julia session we can do
    #prob_loaded = deserialize("my_compiled_problem.jls")
    #sol = solve(prob_loaded)
=#

println("sol time")
@time sol = solve(prob, Tsit5()) #this only takes 10.5 seconds for a 5x5x5 grid
#end

#3x3x3 first run with @mtkcompile: time = 24.6497834, bytes = 2122804680, gctime = 0.7310769, gcstats = Base.GC_Diff(2122804680, 47856, 0, 40144801, 239, 53930, 731076900, 9, 0), lock_conflicts = 0, compile_time = 22.8066422, recompile_time = 2.6308544)
#3x3x3 second run with @mtkcompile: time = 1.7795824, bytes = 354339151, gctime = 0.1219632, gcstats = Base.GC_Diff(354339151, 13193, 0, 7595171, 48, 29811, 121963200, 2, 0), lock_conflicts = 0, compile_time = 0.0003177, recompile_time = 8.78e-5)

#5x5x5 first run with @mtkcompile: time = 148.9056473, bytes = 6397700043, gctime = 3.1595138, gcstats = Base.GC_Diff(6397700043, 148405, 0, 95280353, 234, 140616, 3159513800, 86, 1), lock_conflicts = 0, compile_time = 124.3522357, recompile_time = 0.1753425)
#5x5x5 second run with @mtkcompile: time = 19.6836342, bytes = 4173726251, gctime = 0.7051157, gcstats = Base.GC_Diff(4173726251, 120541, 0, 48697760, 69, 142215, 705115700, 10, 0), lock_conflicts = 0, compile_time = 0.0004347, recompile_time = 0.0004347)

#3x3x3 first run with structural_simplify: time = 17.7417201, bytes = 1749498877, gctime = 0.5585368, gcstats = Base.GC_Diff(1749498877, 40755, 0, 33566670, 197, 42405, 558536800, 8, 0), lock_conflicts = 0, compile_time = 15.5328464, recompile_time = 0.0240274)
#3x3x3 second run with structural_simplify: time = 2.183446, bytes = 354580897, gctime = 0.108073, gcstats = Base.GC_Diff(354580897, 13191, 0, 7612215, 48, 15720, 108073000, 2, 0), lock_conflicts = 0, compile_time = 0.0001742, recompile_time = 0.0001742)

#5x5x5 first run with structural_simplify: time = 133.4732483, bytes = 6389232449, gctime = 2.1426962, gcstats = Base.GC_Diff(6389232449, 148343, 0, 95118570, 233, 148297, 2142696200, 93, 1), lock_conflicts = 0, compile_time = 110.5543221, recompile_time = 0.0008988)
#5x5x5 second run with structural_simplify: time = 20.9268831, bytes = 4174243188, gctime = 1.0864351, gcstats = Base.GC_Diff(4174243188, 120540, 0, 48758377, 69, 132999, 1086435100, 80, 0), lock_conflicts = 0, compile_time = 0.0005768, recompile_time = 0.0005768)

#10x10x1 first run with structural_simplify: time = 66.4354146, bytes = 3994300289, gctime = 1.68664, gcstats = Base.GC_Diff(3994300289, 94479, 0, 65909926, 227, 95745, 1686640000, 46, 1), lock_conflicts = 0, compile_time = 54.0502099, recompile_time = 0.0008226)
#the above took 25 seconds in the 2D FEM with mtkcomponents file (10x10 grid) on first run
#the above took  44.3 seconds in the 3D FEM with mtkcomponents file (10x10x1 grid) on first run. complete output: 44.315889 seconds (58.30 M allocations: 3.390 GiB, 4.98% gc time, 78.12% compilation time: <1% of which was recompilation)
#10x10x1 second run with structural_simplify: time = 14.7961221, bytes = 2180224051, gctime = 0.8638261, gcstats = Base.GC_Diff(2180224051, 66885, 0, 29139042, 64, 75347, 863826100, 29, 0), lock_conflicts = 0, compile_time = 0.0006752, recompile_time = 0.0006752)


#4x4x4 detailed
#=
structural_simplify time
    1.881339 seconds (6.00 M allocations: 405.779 MiB, 9.94% gc time, 0.02% compilation time: 100% of which was recompilation)
prob time
    41.985048 seconds (28.07 M allocations: 1.588 GiB, 1.74% gc time, 80.81% compilation time)
sol time
    8.857977 seconds (20.60 M allocations: 1.065 GiB, 4.18% gc time, 99.98% compilation time)
=#

#4x4x4 detailed second run
#=
structural_simplify time
    3.136403 seconds (5.96 M allocations: 404.536 MiB, 43.32% gc time, 0.01% compilation time: 100% of which was recompilation)
prob = ODEProblem() time
    6.153364 seconds (15.53 M allocations: 1.091 GiB, 5.71% gc time)
sol time
    0.001153 seconds (16.58 k allocations: 566.454 KiB)
=#