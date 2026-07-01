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
        
        for cell in CellIterator(grid)
            dist_between_points = zeros(0.0, 7)
            cell_id = cellid(cell)
            start_coord = getcoordinates(grid, cell_id)[1] 
            #this isn't exactly ideal because we're just measuring the distance bewteen one arbitrary start point 
            #and all the other points but it works for now 
            Ferrite.reinit!(cv, cell)
            for j in 2:getnquadpoints(cv)
                point_coord = getcoordinates(grid, cell_id)[j]
                dist_between_points = norm(start_coord - point_coord)
            end
            
            aspect_ratio = maximum(dist_between_points) / minimum(dist_between_points)
            if aspect_ratio >= max_acceptable_aspect_ratio  
                println("cell $(cell_id)'s aspect ratio of ($aspect_ratio) is greater than max acceptable aspect ratio ($max_acceptable_aspect_ratio)")
            end
        end
    end
end

check_mesh(grid, 500) #For some reason this takes forever
