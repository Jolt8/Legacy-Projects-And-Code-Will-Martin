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

#=
function check_mesh(grid, max_acceptable_aspect_ratio)
    cell_geometries = get_cell_geometries(grid)

    ref_shape = getrefshape(grid.cells[1])
    cell_qr = QuadratureRule{ref_shape}(2)
    poly_interp = Lagrange{ref_shape, 1}() 
    cv = CellValues(cell_qr, poly_interp)

    for (i, conn) in enumerate(cell_geometries)
        vol = conn.volume
        areas = conn.face_areas
        centroid_coords = conn.centroid_coords

        if vol < 0 
            println("vol at $i ($vol) is less than zero")
        elseif vol == 0 
            println("vol at $i ($vol) is zero")
        end
        #TODO: check if aspect ratio is too small
        
        for cell in CellIterator(grid)
            dist_between_points = Float64[]
            cell_id = cellid(cell)
            start_coord = getcoordinates(grid, cell_id)[1]
            Ferrite.reinit!(cv, cell)
            for j in 2:getnquadpoints(cv)
                point_coord = getcoordinates(grid, cell_id)[j]
                push!(dist_between_points, norm(start_coord - point_coord))
            end
            
            aspect_ratio = maximum(dist_between_points) / minimum(dist_between_points)
            if aspect_ratio >= max_acceptable_aspect_ratio  
                println("cell $(cell_id)'s aspect ratio of ($aspect_ratio) is greater than max acceptable aspect ratio ($max_acceptable_aspect_ratio)")
            end
        end
    end
end
=#
#check_mesh(grid, 500) #For some reason this takes forever

#Ferrite.cellcenter(RefHexahedron, Float64) #currently this doesn't have an implementation for pyramids :(


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

function setup_matrix_problem(grid, cell_geometries, cell_sets, facet_sets)
    n_cells = getncells(grid)

    p = setup_manual_problem_all_func(grid, cell_geometries, cell_sets, facet_sets)

    # 2. Create storage for the conductivity matrix K
    I = Int[]
    J_idx = Int[] 
    V = Float64[]
    sizehint!(I, n_cells * 7) # Approx. 7 non-zero entries per row (1 diagonal + 6 neighbors)
    sizehint!(J_idx, n_cells * 7)
    sizehint!(V, n_cells * 7)

    # 3. Populate K from the connections
    for conn in p.connections
        i = conn.cell_idx_a
        j = conn.cell_idx_b
        G = conn.transmissibility # Transmissibility (k*A/d)

        # The term G * (u_j - u_i) contributes to the i-th row.
        # This means K[i,i] += G and K[i,j] -= G.
        # Similarly for the j-th row.

        # Add G to diagonal entries K[i,i] and K[j,j]
        push!(I, i); push!(J_idx, i); push!(V, G)
        push!(I, j); push!(J_idx, j); push!(V, G)

        # Add -G to off-diagonal entries K[i,j] and K[j,i]
        push!(I, i); push!(J_idx, j); push!(V, -G)
        push!(I, j); push!(J_idx, i); push!(V, -G)
    end

    # Assemble the sparse conductivity matrix K
    K = sparse(I, J_idx, V, n_cells, n_cells)

    # 4. Create the source vector `b` (in Watts)
    b = zeros(Float64, n_cells)
    
    # Volumetric sources from NeuBC on cells
    for i in 1:n_cells
        if typeof(p.cell_type_map[i]) == NeuBC
            b[i] += p.cell_type_map[i].flux
        end
        # Note: We can add RobBC handling here later if needed
    end

    # Surface flux from FacetNeuBC on facets
    
    for cell in CellIterator(grid)
        i = cellid(cell)
        for j in 1:nfacets(cell)
            if typeof(p.facet_type_map[i, j]) == FacetNeuBC
                # The flux is per unit area, but your setup seems to treat it as total flux.
                # If it's flux per area, you'd multiply by face area. Assuming total flux for now.
                b[i] += p.facet_type_map[i, j].flux
            end
        end
    end

    # 5. Handle Dirichlet (fixed temperature) boundary conditions
    C_inv_diag = 1.0 ./ p.capacities
    for i in 1:n_cells
        if typeof(p.cell_type_map[i]) == DirBC
            # For a fixed cell, du/dt must be 0.
            # We achieve this by making its corresponding row in the final J and f_source zero.
            C_inv_diag[i] = 0.0 # This zeros the row in -C⁻¹*K
            b[i] = 0.0          # This zeros the source term for this cell
        end
    end
    
    # 6. Create the final system matrix and source vector for the ODE
    # du/dt = -C⁻¹*K*u + C⁻¹*b
    C_inv = spdiagm(0 => C_inv_diag)
    
    J_system = C_inv * -K
    f_source = C_inv * b

    return J_system, f_source, p
end

function FVM_matrix_f!(du, u, p, t)
    J, f_source = p # Unpack the pre-computed matrix and source vector

    # This is the in-place, non-allocating version of: du = J * u + f_source
    mul!(du, J, u)
    du .+= f_source
end

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
copper = Material(401.0u"W/(m*K)", 8.96u"g/cm^3", 0.385u"J/(g*K)") 
steel = Material(45.0u"W/(m*K)", 7.85u"g/cm^3", 0.446u"J/(g*K)")

left_set = CellSet("left", heat_volume_transmissibility, heat_volume_capacity, copper, NeuBC(500.0))
default_set = CellSet("default", heat_volume_transmissibility, heat_volume_capacity, copper, FreeBC())

cell_sets = [left_set, default_set]

#Face BCS Definitions
top = ExclusiveTopology(grid)
addboundaryfacetset!(grid, top, "flux_in", x -> x[1] ≈ 0.0)
getfacetset(grid, "flux_in")
flux_in_set = FacetSet("flux_in", FacetNeuBC(5000.0))

facet_sets = [flux_in_set]

@time J, f_source, p_matrix = setup_matrix_problem(grid, cell_geometries, cell_sets, facet_sets)

ode_params = (J, f_source)

n_cells = getncells(grid)

u0 = fill(293.15, n_cells)

for cell in CellIterator(grid)
    cell_idx = cellid(cell)
    if cell_idx in collect(getcellset(grid, "left"))
        u0[cell_idx] = 500.13
    else
        u0[cell_idx] = 273.13
    end
end

tspan = (0.0, 100.0)

println("Creating matrix-based ODE problem...")

function update_jacobian!(J_out, u, p, t)
    J_precalc, _ = p
    # In a linear constant system, J doesn't change. 
    # We copy it to ensure the solver uses the correct values.
    J_out .= J_precalc 
    return nothing
end

func = ODEFunction(FVM_matrix_f!; jac=update_jacobian!, jac_prototype=J)

prob_matrix = ODEProblem(func, u0, tspan, ode_params)

println("Solving with matrix function (using a stiff solver)...")
desired_amount_of_u = 100
@time sol = solve(prob_matrix, Rodas5(), saveat=(tspan[end]/desired_amount_of_u))

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


ref_shape = getrefshape(grid.cells[1])
poly_interp = Lagrange{ref_shape, 1}() #1 = linear elements 2 = quadratic/curved edges
cell_qr = QuadratureRule{ref_shape}(2) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation
facet_qr = FacetQuadratureRule{ref_shape}(2)

fv = FacetValues(facet_qr, poly_interp)

dh = DofHandler(grid)

iv = InterfaceValues(facet_qr, poly_interp)
#for use later when we require normals
#=
for cell in CellIterator(grid)
    cell_id = cellid(cell)
    for fc in FacetIterator(grid, getfacetset(grid, "flux_in"))
        Ferrite.reinit!(fv, fc)
        for q_point in 1:getnquadpoints(iv)
            normal = getnormal(fv, q_point)
            println(normal)
        end
    end
end=#