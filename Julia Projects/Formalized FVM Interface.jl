using Ferrite
using SparseArrays

struct CellGeometry
    volume::Float64
    face_areas::Vector{Float64}
    centroid_coords::Vector{Float64}
end

function FVM_cell_geometries(grid, poly_interp, cell_quadrature_rule, facet_quadrature_rule)
    cell_qr = cell_quadrature_rule
    cell_values = CellValues(cell_qr, poly_interp)
    facet_qr = facet_quadrature_rule
    facet_values = FacetValues(facet_qr, poly_interp)

    cell_geometries = []

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
        x_avg = sum([c[1] for c in cell_coords]) / length(cell_coords)
        y_avg = sum([c[2] for c in cell_coords]) / length(cell_coords)
        z_avg = sum([c[3] for c in cell_coords]) / length(cell_coords)

        centroid_coords = [x_avg, y_avg, z_avg]

        push!(cell_geometries, CellGeometry(vol_i, areas, centroid_coords)) #formatting it this way because it's kindof 3D, 2D, 1D although a different format might be better
    end
    return cell_geometries
end

#START of generating grid
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

grid = generate_grid(Hexahedron, (2, 2, 2), left, right)
#END of generating grid

cell_type = getrefshape(eltype(grid.cells))
#eltype(grid.cells) #RefHexahedron #TODO: this should maybe be integrated into the function itself and just assume linear and 2 for integration points

poly_interp = Lagrange{base_shape, 1}() #1 = linear elements 2 = quadratic/curved edges

cell_qr = QuadratureRule{base_shape}(2) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation

facet_qr = FacetQuadratureRule{base_shape}(2) 

FVM_grid_dimensions(grid, poly_interp, cell_qr, facet_qr)

#END of main components for face area and cell volume


#START of FVM neighbor finder
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

left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

grid = generate_grid(Hexahedron, (2, 2, 1), left, right)

test = FVM_unique_connections(grid)

println(test) #a 2x2x1 Hexahedron gives 4 unique connections so this definitely works #a 2x2x2 gives 12 elements, which is also right!
#END of FVM neighbor finder


#Some structs that may be useful in the future:
struct VertexConnection
    vertex1_idx::Int
    vertex2_idx::Int
    coords1::Ferrite.Vec{2, Float64}  #or Vec{3, Float64} for 3D
    coords2::Ferrite.Vec{2, Float64}
    distance::Float64
    cross_section_area::Float64  #For 2D: element thickness × perpendicular width
end

struct FacetConnection
    cell1_idx::Int
    cell2_idx::Int
    facet_area::Float64
    facet_normal::Ferrite.Vec{2, Float64}  # Points from cell1 to cell2
    distance_between_centroids::Float64
end

#mtk component matching
struct FEMNodeMapping{T}
    vertex_idx::Int
    component::T  #ex. HeatCapacitor
end

struct FEMEdgeMapping{T}
    connection::VertexConnection
    component::T  #ex. ThermalConductor
end

struct FVMCellMapping{T}
    cell_idx::Int
    component::T
    volume::Float64
end

struct FVMFacetMapping{T}
    connection::FacetConnection
    component::T  #ex. some flux component
end