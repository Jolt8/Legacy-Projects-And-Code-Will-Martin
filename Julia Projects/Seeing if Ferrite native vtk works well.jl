using Ferrite
using WriteVTK

left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

grid_dimensions = (2, 2, 2)

grid = generate_grid(Hexahedron, grid_dimensions, left, right)

VTKGridFile("C://Users//wille//Desktop//Legacy-Projects-And-Code-Will-Martin//Julia Projects//test_solution", grid) do vtk
    write_cell_data(vtk, 1.0, "T")
end
