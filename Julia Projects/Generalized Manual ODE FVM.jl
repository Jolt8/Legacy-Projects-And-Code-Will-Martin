using Ferrite
using DifferentialEquations
using LinearAlgebra #for norm()
using WriteVTK
using Dates
using Unitful

#START OF FERRITE SECTION
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

#grid_dimensions = (4, 4, 4)
grid_dimensions = (100, 10, 10) 
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

addcellset!(grid, "left", x -> x[1] <= left[1] + right[1] / 2)
getcellset(grid, "left")

struct CellGeometry
    volume::typeof(1.0u"m^3")
    face_areas::Vector{typeof(1.0u"m^2")}
    centroid_coords::Vector{3, Float64}
end

function get_cell_geometries(grid)
    poly_interp = Lagrange{getrefshape(grid.cells[1]), 1}() #1 = linear elements 2 = quadratic/curved edges
    cell_qr = QuadratureRule{getrefshape(grid.cells[1])}(2) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation
    facet_qr = FacetQuadratureRule{getrefshape(grid.cells[1])}(2)
    #we should maybe move this to the for cell loop to make sure if it gets a different ref shape from a cell it doesn't break 
    #although we really shouldn't if that takes too much computation
    #not even sure if hex and tetra meshes are ever combined 
    cell_values = CellValues(cell_qr, poly_interp)
    facet_values = FacetValues(facet_qr, poly_interp)

    cell_geometries = Vector{CellGeometry}()

    for cell in CellIterator(grid)
        Ferrite.reinit!(cell_values, cell) #the Ferrite. is required here because reint! is defined by many other packages
        vol = sum(getdetJdV(cell_values, qp) for qp in 1:getnquadpoints(cell_values)) * 1.0u"m^3"
        
        cell_coords = getcoordinates(cell)
        areas = []
        
        for face_idx in 1:nfacets(cell)
            Ferrite.reinit!(facet_values, cell_coords, face_idx) 
            push!(areas, sum(getdetJdV(facet_values, qp) for qp in 1:getnquadpoints(facet_values)) * 1.0u"m^2")
        end

        centroid_vec = sum(cell_coords) / length(cell_coords)
        centroid_coords = [centroid_vec[1], centroid_vec[2], centroid_vec[3]]
        
        push!(cell_geometries, CellGeometry(vol, areas, centroid_coords))
    end
    return cell_geometries
end

cell_geometries = get_cell_geometries(grid)
# END OF FERRITE SECTION

#Ideally I'd like to make writing these kind of things as similar to mtk components as possible 

#A more generalized version
struct Connection
    cell_idx_a::Int
    cell_idx_b::Int
    transmissibility::typeof(1.0u"J/K") #not sure what conductance would be generalized to here
    #perhaps we should also store distance from cell centers here
end

# A struct to hold all parameters for the ODE problem
struct HeatProblemParameters
    capacities::Vector{typeof(1.0u"J/K")}
    connections::Vector{Connection}
    heat_sources::Vector{Tuple{Int, typeof(1.0u"W")}} #(cell_index, Q_flow)
end

struct CellGroupDefinition
    ferrite_name::AbstractString
    type::AbstractString
    parameters::Vector{Dict{Int, Vector{Any}}}
    #parameters::Vector{Tuple{Int, Vector{Any}}} #not sure if a dict, vector, or tuple is better here
end
#not sure if the bottom of top is better
struct CellGroupDefinition
    ferrite_name::AbstractString
    type::AbstractString
    parameters::Vector{Any}
end 



#we could either make types for Dirichet, Neumann, Robin and normal cell or
#we could make cell types for bcs and non_bcs
CellGroupDefinition("left", "bcs", [Dict(1 => [273.13u"K"])])
#CellGroupDefinition("left", "bcs", [(1, [273.13u"K"])]) 
CellGroupDefinition("unnamed", "non_bcs", [Dict(1 => [273.13u"K"])]) #unnamed/internal
#or
CellGroupDefinition("left", "Dirichet", [Dict(1 => [273.13u"K"])])
CellGroupDefinition("unnamed", "non_bcs", [Dict(1 => [273.13u"K"])])


struct CellGroupDefinition
    ferrite_name::AbstractString
    type::AbstractString
    parameters::Vector{Any}
end


struct ProblemParameters
    capacities::Vector{typeof(1.0u"J/K")} #capacitance/flux resistance can't come up with a good general name
    connections::Vector{Connection}
    cell_types::Vector{CellType} 
end


# 2. CREATE A SETUP FUNCTION

function setup_manual_problem(grid, cell_geometries) 
    #k_thermal, rho and cp should be replaced by a struct that takes in different parameters for each cell type/bcs
    n_cells = getncells(grid)
    vol_heat_cap = rho * cp

    #calculate heat capacity in J/K for each cell
    capacities = [geo.volume * vol_heat_cap for geo in cell_geometries]

    connections = Connection[]
    top = ExclusiveTopology(grid)

    for i in 1:n_cells
        for j in 1:nfacets(grid.cells[i])
            neighbor_info = top.face_face_neighbor[i, j]
            if !isempty(neighbor_info)
                neighbor_idx = collect(neighbor_info[1].idx)[1] #the collect()[1] is required here because neighbor_info[1].idx returns a tuple 
                
                #IMPORTANT: To avoid double-counting connections (e.g., 1->2 and 2->1),
                # we only process the connection if our cell index is smaller.
                if i < neighbor_idx
                    face_area = cell_geometries[i].face_areas[j]
                    dist = norm(cell_geometries[i].centroid_coords - cell_geometries[neighbor_idx].centroid_coords) * 1.0u"m"
                    G_calc = k_thermal * face_area / dist
                    
                    push!(connections, Connection(i, neighbor_idx, G_calc))
                end
            end
        end
    end
    
    #define bcs
    for flux_bcs in flux_bcs_
        flux_bcs_node_indicies = collect(getcellset(flux_bcs.ferrite_name))

        #left_node_indices = collect(getcellset(grid, "left"))
        #heat_sources = [(idx, 5000.0u"W") for idx in left_node_indices]
    end 
    return ProblemParameters(capacities, connections)
end


# 3. WRITE THE ODE FUNCTION `f!`

function FVM_iter_f!(du, u, p::HeatProblemParameters, t)
    # `du` is the vector of dT/dt for each cell
    # `u` is the vector of current temperatures T for each cell
    # `p` is our struct holding all the system parameters
    # `t` is the current time (not used in this problem, but required by the solver)

    du .= 0.0 * Unitful.unit(u) / 1.0u"s"

    #Calculate heat flow from internal connections
    #Q_flow = G * (T_a - T_b)
    for conn in p.connections
        i = conn.cell_idx_a
        j = conn.cell_idx_b
        G = conn.conductance

        temp_diff = u[i] - u[j]
        flux = G * temp_diff

        #the flux leaves cell i and enters cell j
        du[i] -= flux / p.capacities[i]
        du[j] += flux / p.capacities[j]
    end

    #add heat flow from boundary conditions (heat sources)
    #=for (idx, Q_flow) in p.heat_sources
        du[idx] += Q_flow / p.capacities[idx]
    end=#

    #Step D: Convert net heat flow (Q_net, which is currently in `du`)to rate of temperature change (dT/dt) using the formula:
    #C * dT/dt = Q_net  =>  dT/dt = Q_net / C
    #du ./= p.capacities 
    #this is actually a limitation of unitful because we can't do du ./= p.capacities because if we don't do that while calculating du, we get a unit error
    #we could get around this with a pre_du vector that's divided to yield du afterwards but that might be more computationally expensive
end


# 4. SETUP AND SOLVE

# Material properties for copper
k_thermal = 400u"W/(m*K)"
rho = 9000.0u"kg/m^3"
cp = 0.385u"J/(g*K)"

println("Setting up manual problem...")
@time p = setup_manual_problem(grid, cell_geometries)

#setup
n_cells = getncells(grid)

#n_left_bcs = length(getcellset(grid, "left"))

u0 = fill(293.15u"K", n_cells)

left_bcs_cell_idxs = collect(getcellset(grid, "left"))
left_bcs_cell_temp = 500.13u"K"

regular_cell_idxs = "something"
regular_cell_temp = 273.13u"K"

for cell in CellIterator(grid)
    cell_idx = cellid(cell)
    if cell_idx in left_bcs_cell_idxs
        u0[cell_idx] = left_bcs_cell_temp
    else
        u0[cell_idx] = regular_cell_temp
    end
end

u0

tspan = (0.0u"s", 100u"s")

println("prob = ODEProblem() time (manual)")

@time prob = ODEProblem(heat_transfer_f!, u0, tspan, p)

println("sol time (manual)")
@time sol = solve(prob, Tsit5(), dtmin=0.001) #I wonder if there's a way to automatically adjust the time span to get around x iterations
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
record_sol = false

if record_sol == true
    date_and_time = Dates.format(now(), "yyyy-mm-dd I.MM.SS p")
    date_and_time = Dates.format(now(), "I.MM.SS p")

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
        state_at_t = ustrip.(sol.u[step])
        t_stripped = ustrip(t)
        t_rounded = ustrip(round(typeof(1.0u"s"), t, digits=4))
        VTKGridFile(step_filename * " $step" * " at $t_rounded.vtu", grid) do vtk
            write_cell_data(vtk, state_at_t, "T")
            pvd[t_stripped] = vtk
        end
    end
    vtk_save(pvd)
end
