using Ferrite
using WriteVTK

#check out polyhedron cells: https://juliavtk.github.io/WriteVTK.jl/dev/grids/unstructured/#Unstructured-grid
#here's one for multiple cells: https://github.com/JuliaVTK/WriteVTK.jl/blob/master/test/polyhedron_cube.jl 

left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

grid_dimensions = (2, 2, 2)

grid = generate_grid(Hexahedron, grid_dimensions, left, right)

points = Float32[]
cells = VTKPolyhedron[]
for cell in CellIterator(grid)
    cell_idx = cellid(cell)
    
    node_coords = getcoordinates(grid, cell_idx)
    append!(points, collect(Iterators.flatten(node_coords)))

    nodes_of_cell = collect(grid.cells[cell_idx].nodes)
    
    nodes_of_cell_faces = Ferrite.faces(grid.cells[cell_idx])

    println(nodes_of_cell_faces)

    push!(cells, VTKPolyhedron(nodes_of_cell, nodes_of_cell_faces))
end
points = reshape(points, 3, :)

points
cells

vtk_grid("polyhedron_cube", points, cells; compress = false) do vtk
    vtk["T"] = vec(sum(x -> abs2(x + 1), points; dims = 1))
end

vpoints = Float32[
    -1, -1, -1,
     1, -1, -1,
     1,  1, -1,
    -1,  1, -1,
    -1, -1,  1,
     1, -1,  1,
     1,  1,  1,
    -1,  1,  1,
]

append!(vpoints, vpoints .+ 2)  # append points for second cube
points = reshape(vpoints, 3, :)

cells = [
    # Cube 1
    VTKPolyhedron(
        1:8,           # connectivity vector
        (1, 4, 3, 2),  # face 1
        (1, 5, 8, 4),  # face 2
        (5, 6, 7, 8),  # etc...
        (6, 2, 3, 7),
        (1, 2, 6, 5),
        (3, 4, 8, 7),
    ),
    # Cube 2
    VTKPolyhedron(
        8 .+ (1:8),         # connectivity vector
        8 .+ (1, 4, 3, 2),  # face 1
        8 .+ (1, 5, 8, 4),  # face 2
        8 .+ (5, 6, 7, 8),  # etc...
        8 .+ (6, 2, 3, 7),
        8 .+ (1, 2, 6, 5),
        8 .+ (3, 4, 8, 7),
    ),
]
points

@time vtk_grid("polyhedron_cube", points, cells; compress = false) do vtk
    vtk["temperature"] = vec(sum(x -> abs2(x + 1), points; dims = 1))
    vtk["pressure"] = [42.0, 84.0]  # one per cube
end