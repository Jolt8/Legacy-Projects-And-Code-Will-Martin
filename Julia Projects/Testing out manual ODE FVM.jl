using Ferrite
using DifferentialEquations
using LinearAlgebra # For norm()

# Keep the Ferrite part exactly as it is to generate the grid
# and cell geometries. I'll copy the relevant parts for a self-contained example.
# START OF FERRITE SECTION
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

#grid_dimensions = (4, 4, 4)
grid_dimensions = (100, 100, 100) 
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

length_to_node_ratio = right[1] / collect(grid_dimensions)[1]

addcellset!(grid, "left", x -> x[1] <= left[1] + length_to_node_ratio)
addcellset!(grid, "right", (x) -> x[1] >= (right[1] - 0.0000001) - (length_to_node_ratio))

struct CellGeometry
    volume::Float64
    face_areas::Vector{Float64}
    centroid_coords::Vector{Float64}
end

# This function is simplified to just get what we need for the manual setup
function get_cell_geometries(grid)
    poly_interp = Lagrange{RefHexahedron, 1}()
    cell_qr = QuadratureRule{RefHexahedron}(2)
    facet_qr = FacetQuadratureRule{RefHexahedron}(2)
    cell_values = CellValues(cell_qr, poly_interp)
    facet_values = FacetValues(facet_qr, poly_interp)

    cell_geometries = Vector{CellGeometry}()

    for cell in CellIterator(grid)
        Ferrite.reinit!(cell_values, cell)
        vol = sum(getdetJdV(cell_values, qp) for qp in 1:getnquadpoints(cell_values))
        
        cell_coords = getcoordinates(cell)
        areas = [begin
                    Ferrite.reinit!(facet_values, cell_coords, face_idx)
                    sum(getdetJdV(facet_values, qp) for qp in 1:getnquadpoints(facet_values))
                end for face_idx in 1:nfacets(cell)]

        centroid_vec = sum(cell_coords) / length(cell_coords)
        centroid_coords = [centroid_vec[1], centroid_vec[2], centroid_vec[3]]

        push!(cell_geometries, CellGeometry(vol, areas, centroid_coords))
    end
    return cell_geometries
end

cell_geometries = get_cell_geometries(grid)
# END OF FERRITE SECTION


# 1. DEFINE DATA STRUCTURES

# A struct to hold info about a connection between two cells
struct ThermalConnection
    cell_idx_a::Int
    cell_idx_b::Int
    conductance_G::Float64
end

# A struct to hold all parameters for the ODE problem
struct HeatProblemParameters
    capacities::Vector{Float64}
    connections::Vector{ThermalConnection}
    heat_sources::Vector{Tuple{Int, Float64}} # (cell_index, Q_flow)
end


# 2. CREATE A SETUP FUNCTION

function setup_manual_problem(grid, cell_geometries, k_thermal, rho, cp)
    n_cells = getncells(grid)
    vol_heat_cap = rho * cp

    # Calculate heat capacity for each cell
    capacities = [geo.volume * vol_heat_cap for geo in cell_geometries]

    # Build the list of unique connections and their conductances
    connections = ThermalConnection[]
    top = ExclusiveTopology(grid) # Helps find neighbors efficiently

    for i in 1:n_cells
        for j in 1:nfacets(grid.cells[i])
            neighbor_info = top.face_face_neighbor[i, j]
            if !isempty(neighbor_info)
                neighbor_idx = collect(neighbor_info[1].idx)[1]
                
                # IMPORTANT: To avoid double-counting connections (e.g., 1->2 and 2->1),
                # we only process the connection if our cell index is smaller.
                if i < neighbor_idx
                    face_area = cell_geometries[i].face_areas[j]
                    dist = norm(cell_geometries[i].centroid_coords - cell_geometries[neighbor_idx].centroid_coords)
                    G_calc = k_thermal * face_area / dist
                    
                    push!(connections, ThermalConnection(i, neighbor_idx, G_calc))
                end
            end
        end
    end
    
    # Define boundary conditions
    left_node_indices = collect(getcellset(grid, "left"))
    heat_sources = [(idx, 5000.0) for idx in left_node_indices]

    # Package everything into our parameter struct
    return HeatProblemParameters(capacities, connections, heat_sources)
end


# 3. WRITE THE ODE FUNCTION `f!`

function heat_transfer_f!(du, u, p::HeatProblemParameters, t)
    # `du` is the vector of dT/dt for each cell
    # `u` is the vector of current temperatures T for each cell
    # `p` is our struct holding all the system parameters
    # `t` is the current time (not used in this problem, but required by the solver)

    # Step A: Initialize derivatives to zero. This represents starting with no heat flow.
    du .= 0.0

    # Step B: Calculate heat flow from internal connections
    # Q_flow = G * (T_a - T_b)
    for conn in p.connections
        i = conn.cell_idx_a
        j = conn.cell_idx_b
        G = conn.conductance_G

        # Temperature difference drives the flow
        temp_diff = u[i] - u[j]
        flux = G * temp_diff

        # The flux leaves cell `i` and enters cell `j`
        du[i] -= flux
        du[j] += flux
    end

    # Step C: Add heat flow from boundary conditions (heat sources)
    for (idx, Q_flow) in p.heat_sources
        du[idx] += Q_flow
    end

    # Step D: Convert net heat flow (Q_net, which is currently in `du`)
    # to rate of temperature change (dT/dt) using the formula:
    # C * dT/dt = Q_net  =>  dT/dt = Q_net / C
    du ./= p.capacities
end


# 4. SETUP AND SOLVE

# Material properties
k_thermal = 200.0
rho = 2700.0
cp = 900.0

println("Setting up manual problem...")
@time p = setup_manual_problem(grid, cell_geometries, k_thermal, rho, cp)

# --- Simulation setup ---
n_cells = getncells(grid)
u0 = fill(293.15, n_cells) # Initial condition: all cells at 20°C
tspan = (0.0, 10.0)

println("prob = ODEProblem() time (manual)")
# This step is now extremely fast because there is no code generation!
@time prob_manual = ODEProblem(heat_transfer_f!, u0, tspan, p)

println("sol time (manual)")
@time sol_manual = solve(prob_manual, Tsit5())

println("\nDone. Notice how the `ODEProblem` creation time is now negligible.")