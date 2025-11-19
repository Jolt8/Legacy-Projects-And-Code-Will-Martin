using Ferrite
using SparseArrays
#START of generating grid
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

grid = generate_grid(Hexahedron, (2, 2, 2), left, right)
#END of generating grid

#START of main components for face area and cell volume
ip = Lagrange{RefHexahedron, 1}() #no idea what the 1 represents here

cell_qr = QuadratureRule{RefHexahedron}(2)
cell_values = CellValues(cell_qr, ip)

facet_qr = FacetQuadratureRule{RefHexahedron}(2) #no idea what the 2 represents here
facet_values = FacetValues(facet_qr, ip)

#cell_coords = getcoordinates(grid, 1) #this returns a cube with 1.0 long sides 

cell_dimensions = []

"""
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
cell_dimensions
"""

function FVM_grid_dimensions(grid, poly_interp, cell_quadrature_rule, facet_quadrature_rule)
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
        areas = []
        cell_coords = getcoordinates(cell)
        
        for facet_index in 1:nfacets(cell)
            reinit!(facet_values, cell_coords, facet_index)
            
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

        push!(cell_dimensions, [vol_i, [areas], centroid_coords]) #formatting it this way because it's kindof 3D, 2D, 1D although a different format might be better
    end
    return cell_dimensions
end

poly_interp = Lagrange{RefHexahedron, 1}() #1 = linear elements 2 = quadratic/curved edges

cell_qr = QuadratureRule{RefHexahedron}(2) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation

facet_qr = FacetQuadratureRule{RefHexahedron}(2) 

FVM_grid_dimensions(grid, poly_interp, cell_qr, facet_qr)
#END of main components for face area and cell volume

#START of FVM neighbor finder
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

grid = generate_grid(Hexahedron, (2, 2, 2), left, right)

top = ExclusiveTopology(grid)

getncells(grid)

top.cell_neighbor

top.face_face_neighbor

top.face_face_neighbor[1, 5]
top.face_face_neighbor[1, 3]
top.face_face_neighbor[2, end] 
#If cell 1 had connection with cell 2 at (cell 1, face 3) and (cell 2, face 5)
#OH, I get how this is formatted, using the example above, 1 is the first index of the matrix for cell 1, 
#the 3rd index of the vector within the array is the index of the corresponding face idx for cell 1,
#that vector at the index [1, 3] of the matrix contains the struct FaceIndex that contains (2, 5)
#this means that cell 2 has a connection with cell 1 at cell 2's 5th side 

#this is confirmed by the fact that at index 5 (5th side) of the 2nd index (cell 2) in the matrix contains (FaceIndex(1, 3) or a connection with face 3 of cell 1

unique_connections = Set{Tuple{Tuple{Number, Number}, Tuple{Number, Number}}}() #we should use a set but I don't know how to set it up right

top.face_face_neighbor[1, 3]

for cell in CellIterator(grid)
    i = cellid(cell) 
    ff_neighbor_matrix = top.face_face_neighbor
    for j in 1:nfacets(cell)
        if isempty(ff_neighbor_matrix[i, j])
            nothing
        else
            #[i, j] = [1, 3] #[cell, face]
            #ff_neighbor_matrix[i, j] = [FacetIndex(2, 5)]
            #ff_neighbor_matrix[i, j][1].idx = (2, 5)
            #collect(ff_neighbor_matrix[i, j][1].idx) = [2, 5] #[cell, face]
            #neighbor_cell_and_face = collect(ff_neighbor_matrix[i, j][1].idx)
            neighbor_cell_idx = collect(ff_neighbor_matrix[i, j][1].idx)[1] #I don't know which of these is clearer (top vs bottom two)
            connecting_face_idx = collect(ff_neighbor_matrix[i, j][1].idx)[2]

            connection = minmax((i, j), (neighbor_cell_idx, connecting_face_idx))
            push!(unique_connections, connection)
        end
    end
end
unique_connections

#cleaned up version
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

test = FVM_unique_connections(grid)

println(test)#a 2x2x1 Hexahedron gives 4 unique connections so this definitely works #a 2x2x2 gives 12 elements, which is also right!


CellIndex(1)

facetskeleton(top, grid)

for cell in CellIterator(grid)
    cell_idx = CellIndex(cellid(cell))
    neighborhood = getneighborhood(top, grid, cell_idx, false)
    println(neighborhood)
    for i in 1:nfacets(cell)
        test = getneighborhood(top, grid, FacetIndex((cellid(cell), i)))
        #println(test)
    end
    for facet_idx in facetskeleton(top, grid) #this is formatted FacetIndex(cell_idx, side_number) #ex FacetIndex(1, 5)
    end
end

#END of FVM neighbor finder




#Some structs that may be useful in the future:
struct VertexConnection
    vertex1_idx::Int
    vertex2_idx::Int
    coords1::Vec{2, Float64}  #or Vec{3, Float64} for 3D
    coords2::Vec{2, Float64}
    distance::Float64
    cross_section_area::Float64  #For 2D: element thickness × perpendicular width
end

struct FacetConnection
    cell1_idx::Int
    cell2_idx::Int
    facet_area::Float64
    facet_normal::Vec{2, Float64}  # Points from cell1 to cell2
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