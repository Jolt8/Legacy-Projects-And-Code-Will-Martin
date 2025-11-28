using Ferrite
using DifferentialEquations
using LinearAlgebra #for norm()
using WriteVTK
using Dates
using Unitful
using SparseArrays

#START OF FERRITE SECTION
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

#grid_dimensions = (4, 4, 4)
grid_dimensions = (100, 10, 10) 
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

addcellset!(grid, "left", x -> x[1] <= left[1] + right[1] / 2)
getcellset(grid, "left")
addcellset!(grid, "default", x -> x[1] >= left[1] + right[1] / 2)
getcellset(grid, "default")

struct CellGeometry
    volume::Float64
    face_areas::Vector{Float64}
    centroid_coords::Vector{Float64}
end

function get_cell_geometries(grid)
    ref_shape = getrefshape(grid.cells[1])
    poly_interp = Lagrange{ref_shape, 1}() #1 = linear elements 2 = quadratic/curved edges
    cell_qr = QuadratureRule{ref_shape}(2) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation
    facet_qr = FacetQuadratureRule{ref_shape}(2)
    #we should maybe move this to the for cell loop to make sure if it gets a different ref shape from a cell it doesn't break 
    #although we really shouldn't if that takes too much computation
    #not even sure if hex and tetra meshes are ever combined 
    cell_values = CellValues(cell_qr, poly_interp)
    facet_values = FacetValues(facet_qr, poly_interp)

    cell_geometries = Vector{CellGeometry}()

    for cell in CellIterator(grid)
        Ferrite.reinit!(cell_values, cell) #the Ferrite. is required here because reint! is defined by many other packages DO NOT OMIT Ferrite.
        vol = sum(getdetJdV(cell_values, qp) for qp in 1:getnquadpoints(cell_values))
        
        cell_coords = getcoordinates(cell)
        areas = []
        
        for face_idx in 1:nfacets(cell)
            Ferrite.reinit!(facet_values, cell_coords, face_idx) 
            push!(areas, sum(getdetJdV(facet_values, qp) for qp in 1:getnquadpoints(facet_values)))
        end

        centroid_vec = sum(cell_coords) / length(cell_coords)
        centroid_coords = [centroid_vec[1], centroid_vec[2], centroid_vec[3]]
        
        push!(cell_geometries, CellGeometry(vol, areas, centroid_coords))
    end
    return cell_geometries
end

cell_geometries = get_cell_geometries(grid)
# END OF FERRITE SECTION

abstract type BC end

struct DirBC <: BC
    fixed::Float64
end

struct NeuBC <: BC
    flux::Float64
end

struct RobBC <: BC
    fixed::Float64
    flux::Float64
end

struct FreeBC <: BC
end

abstract type FacetBC end

struct FacetNeuBC <: FacetBC
    flux::Float64
end

struct FacetRobBC <: FacetBC
    fixed::Float64
    flux::Float64
end

struct FacetFreeBC <: FacetBC
end

#These parameters would require units in their base form (no kg) on struct definition but kg can be entered in later
abstract type CellProperties end
struct Material <: CellProperties
    k::typeof(1.0u"W/(m*K)") 
    rho::typeof(1.0u"g/m^3")
    cp::typeof(1.0u"J/(g*K)")
end

struct Connection
    cell_idx_a::Int
    cell_idx_b::Int
    transmissibility::Float64 # k * area / dist
end

struct ProblemParameters
    capacities::Vector{Float64} #capacitance/flux resistance can't come up with a good general name
    connections::Vector{Connection}
    cell_type_map::Vector{BC} 
    facet_type_map::Matrix{FacetBC}
end

# 2. CREATE A SETUP FUNCTION
function setup_manual_problem_all_func(grid, cell_geometries, cell_sets, facet_sets) #another possible version
    n_cells = getncells(grid)
    top = ExclusiveTopology(grid)

    material_id_map = zeros(Int, n_cells)

    #create a map that maps each 
    #FIXME: we need a way to handle cases where some cells aren't assigned to a set, this currently also happens when doing an uneven length grid like (3, 3, 3)
    for (set_index, cell_set) in enumerate(cell_sets)
        set_indices = getcellset(grid, cell_set.set_name)
        for cell_idx in set_indices
            material_id_map[cell_idx] = set_index
        end
    end

    facet_id_map = zeros(Int, n_cells, 6)
    
    for (set_index, facet_set) in enumerate(facet_sets)
        set_indices = collect(getfacetset(grid, facet_set.set_name))
        for facet_index in set_indices #facet_index is formatted FacetIndex(101, 5) where 101 is the cellid and 5 is the faceid
            cell_id = collect(facet_index.idx)[1]
            face_id = collect(facet_index.idx)[2]
            facet_id_map[cell_id, face_id] = set_index
        end
    end

    facet_type_map = Array{FacetBC}(undef, n_cells, 6)

    cell_type_map = Vector{BC}(undef, n_cells)

    capacities = Vector{Float64}(undef, n_cells)

    connections = Connection[]
    sizehint!(connections, n_cells * 3)

    for i in 1:n_cells
        mat_id = material_id_map[i]
        active_set = cell_sets[mat_id]

        #calculate capacity
        cap_val = active_set.capacity_func(cell_geometries[i].volume, active_set.properties)
        capacities[i] = ustrip(cap_val)
        
        #For Dirichlet bcs
        cell_type_map[i] = active_set.type

        for j in 1:nfacets(grid.cells[i])
            facet_id = facet_id_map[i, j]
            if facet_id != 0
                facet_type_map[i, j] = facet_sets[facet_id].type
            else
                facet_type_map[i, j] = FacetFreeBC()
            end

            neighbor_info = top.face_face_neighbor[i, j]
            if !isempty(neighbor_info)
                neighbor_idx = collect(neighbor_info[1].idx)[1] #the collect()[1] is required here because neighbor_info[1].idx returns a tuple 
                
                #IMPORTANT: To avoid double-counting connections (e.g., 1->2 and 2->1),
                # we only process the connection if our cell index is smaller.
                if i < neighbor_idx
                    face_area = cell_geometries[i].face_areas[j]
                    dist = norm(cell_geometries[i].centroid_coords - cell_geometries[neighbor_idx].centroid_coords)
                    
                    neighbor_mat_id = material_id_map[j]
                    neighbor_set = cell_sets[neighbor_mat_id]

                    transmissibility = active_set.transmissibility_func(face_area, dist, active_set.properties, neighbor_set.properties)

                    push!(connections, Connection(i, neighbor_idx, ustrip(transmissibility)))
                end
            end
        end
    end
    return ProblemParameters(capacities, connections, cell_type_map, facet_type_map)
end

#the one thing I don't like about this workflow is that I can't print inside FVM_iter_f! because of the !
function FVM_iter_f!(du, u, p::ProblemParameters, t)
    du .= 0.0

    #Calculate heat flow from internal connections
    #Q_flow = G * (T_a - T_b)
    
    for conn in p.connections
        i = conn.cell_idx_a
        j = conn.cell_idx_b
        G = conn.transmissibility

        temp_diff = u[i] - u[j]
        flux = G * temp_diff

        #the flux leaves cell i and enters cell j
        if typeof(p.cell_type_map[i]) != DirBC
            du[i] -= flux / p.capacities[i]
        end
        if typeof(p.cell_type_map[j]) != DirBC
            du[j] += flux / p.capacities[j]
        end
        if typeof(p.cell_type_map[i]) == DirBC
            du[i] = 0
        end
        if typeof(p.cell_type_map[j]) == DirBC #don't think this is necessary
            du[j] = 0
        end
    end

    for i in 1:length(u)
        if typeof(p.cell_type_map[i]) == NeuBC
            du[i] += p.cell_type_map[i].flux / p.capacities[i]
        elseif typeof(p.cell_type_map[i]) == RobBC
            #something here
        end
    end
    
    for cell in CellIterator(grid)
        i = cellid(cell)
        for j in 1:nfacets(cell)
            if typeof(p.facet_type_map[i, j]) == FacetNeuBC
                du[i] += p.facet_type_map[i, j].flux / p.capacities[i]
            end
        end
    end
    

    #Step D: Convert net heat flow (Q_net, which is currently in `du`)to rate of temperature change (dT/dt) using the formula:
    #C * dT/dt = Q_net  =>  dT/dt = Q_net / C
    #du ./= p.capacities 
    #this is actually a limitation of unitful because we can't do du ./= p.capacities because if we don't do that while calculating du, we get a unit error
    #we could get around this with a pre_du vector that's divided to yield du afterwards but that might be more computationally expensive
end

top = ExclusiveTopology(grid)

# 4. SETUP AND SOLVE
function heat_volume_transmissibility(face_area, dist, p1, p2)
    average_k = (p1.k + p2.k) / 2
    return (average_k * face_area) / dist
end
function heat_volume_capacity(volume, p)
    return p.rho * p.cp * volume
end
#I really like this function idea
#all transmissibilities would require face_area and dist even if they don't use it
#all capacities would require volume

struct CellSet{T, C, P} 
    set_name::String
    transmissibility_func::T
    capacity_func::C
    properties::P
    type::BC #would be dir, neu, rob
end

struct FacetSet 
    set_name::String
    type::FacetBC #would be dir, neu, rob
end

#Cell Definitions
copper = Material(401.0u"W/(m*K)", 8.96u"g/cm^3", 0.385u"J/(g*K)") #FIXME
steel = Material(45.0u"W/(m*K)", 7.85u"g/cm^3", 0.446u"J/(g*K)") #FIXME

left_set = CellSet("left", heat_volume_transmissibility, heat_volume_capacity, copper, NeuBC(500.0))
default_set = CellSet("default", heat_volume_transmissibility, heat_volume_capacity, copper, FreeBC())

cell_sets = [left_set, default_set]

#Face BCS Definitions
addboundaryfacetset!(grid, top, "flux_in", x -> x[1] ≈ 0.0)
getfacetset(grid, "flux_in")
flux_in_set = FacetSet("flux_in", FacetNeuBC(500.0))

facet_sets = [flux_in_set]

#Get p
@time p = setup_manual_problem_all_func(grid, cell_geometries, cell_sets, facet_sets)

p.facet_type_map

#setup
n_cells = getncells(grid)

#n_left_bcs = length(getcellset(grid, "left"))

u0 = fill(293.15, n_cells)

for cell in CellIterator(grid)
    cell_idx = cellid(cell)
    if cell_idx in collect(getcellset(grid, "left"))
        u0[cell_idx] = 500.13
    else
        u0[cell_idx] = 273.13
    end
end

#=
#TODO:
Almost done! - Generate a sparsity pattern to use Implicit Solvers (KenCarp4) for finer grids.
=#
tspan = (0.0, 1.0) #for some reason I can only get change if tspan is absolutely massive

println("prob = ODEProblem() time (manual)")

@time prob = ODEProblem(FVM_iter_f!, u0, tspan, p)

#Defining a sparse array
#=
ref_shape = getrefshape(grid.cells[1])
ip = Lagrange{ref_shape, 1}()
qr = QuadratureRule{ref_shape}(2)

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh);
K = allocate_matrix(dh)
=#

println("sol time (manual)")
desired_amount_of_u = 100
@time sol = solve(prob, Tsit5(), dtmin=0.001, saveat=(tspan[end]/desired_amount_of_u)) 
#I wonder if there's a way to automatically adjust the time span to get around x iterations

println("\n# of u: ", length(sol.u))#using this to roughly estimate how long writing the vtk file will take 

#for a 50x50x50 grid:
#WITH unitful p = setup_manual_problem(grid, cell_geometries, k_thermal, rho, cp) takes 0.800254 seconds (864.07 k allocations: 302.704 MiB, 17.68% compilation time)
#and sol = solve(prob, Tsit5()) takes 32.991329 seconds (3.29 M allocations: 9.556 GiB, 14.99% gc time, 6.51% compilation time: 100% of which was recompilation)

#WITHOUT unitful p = setup_manual_problem(grid, cell_geometries, k_thermal, rho, cp) takes 1.503231 seconds (862.20 k allocations: 302.624 MiB, 51.50% gc time, 7.00% compilation time)
#and sol = solve(prob, Tsit5()) takes 23.891263 seconds (79.42 k allocations: 9.378 GiB, 9.37% gc time, 0.22% compilation time)

#WITHOUT unitful AND WITHOUT du ./= p.capacities p = setup_manual_problem(grid, cell_geometries, k_thermal, rho, cp) takes 0.963267 seconds (862.20 k allocations: 302.624 MiB, 11.73% gc time, 13.97% compilation time)
#and sol = solve(prob, Tsit5()) takes 25.190752 seconds (72.44 k allocations: 9.385 GiB, 10.10% gc time, 0.16% compilation time)

#replace(basename(@__FILE__),r".jl" => "")

#Making the solution avaliable to paraview
record_sol = true

if record_sol == true
    date_and_time = Dates.format(now(), "I.MM.SS p yyyy-mm-dd")
    #date_and_time = Dates.format(now(), "I.MM.SS p")

    root_dir = "C://Users//wille//Desktop//Julia_cfd_output_files"

    project_name = replace(basename(@__FILE__),r".jl" => "")

    #project_name = "HeatSim "

    sim_folder_name = project_name * " " * date_and_time

    output_dir = joinpath(root_dir, sim_folder_name)

    mkpath(output_dir)

    pvd_filename = joinpath(output_dir, "solution_collection")

    pvd = paraview_collection(pvd_filename)

    step_filename = joinpath(output_dir, "timestep")

    #this step may actually become a significant bottleneck
    #update: This is a very big bottleneck
    for (step, t) in enumerate(sol.t)
        state_at_t = sol.u[step]
        VTKGridFile(step_filename * " $step" * " at $t.vtu", grid) do vtk
            write_cell_data(vtk, state_at_t, "T")
            pvd[t] = vtk
        end
    end
    vtk_save(pvd)
end
